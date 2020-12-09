#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include <cmath>
#include <map>
#include <set>

std::vector<int> split(const std::string& line, char separator = ',')
{
    std::vector<int> ret;
    std::string::size_type start_pos{0}, sep_pos;
    while ((sep_pos = line.find(separator, start_pos)) != std::string::npos) {
        ret.push_back(std::stoi(line.substr(start_pos, sep_pos - start_pos)));
        start_pos = sep_pos + 1;
    }
    ret.push_back(std::stoi(line.substr(start_pos)));
    return ret;
}

struct Task
{
    int id, time, question_id, sla_time;
};

struct Assignment
{
    int task_id, processor_id, start_time;
};

int main(int argc, char **argv)
{
    std::string line;

    std::map<int, int> minTime; // (taskid -> minTime)
    std::set<int> processor_ids;
    std::vector<std::vector<int>> pt_matrix;
    std::ifstream ifs_mat("/tcdata/sic_semi_process_time_matrix.csv");
    std::getline(ifs_mat, line);
    pt_matrix.push_back(split(line));
    while (std::getline(ifs_mat, line)) {
        auto v = split(line);
        processor_ids.insert(v[0]);
        pt_matrix.emplace_back(v);

        for (int i = 1; i < v.size(); i++) {
            if (minTime.find(i) == minTime.end()) {
                minTime[i] = v[i];
            } else {
                minTime[i] = std::min(minTime[i], v[i]);
            }
        }
    }

    std::unordered_set<int> task_ids;
    std::unordered_map<int, Task> tasks;
    std::ifstream ifs_order("/tcdata/sic_semi_work_order.csv");
    while (std::getline(ifs_order, line)) {
        auto v = split(line);
        Task task{v[0], v[1], v[2], v[3]};
        task_ids.insert(task.id);
        tasks[task.id] = task;
    }

    std::unordered_set<int> result_task_ids;
    std::unordered_set<int> result_processor_ids;
    std::unordered_map<int, std::vector<Assignment>> result_tasks;
    std::unordered_map<int, std::vector<std::pair<int, int>>> result_processors;
    std::ifstream ifs_result{argv[1]};
    while (std::getline(ifs_result, line)) {
        auto v = split(line);
        int task_id = v[0];
        int processor_id = v[1];
        int start_time = v[2];
        result_task_ids.insert(task_id);
        result_processor_ids.insert(processor_id);
        result_tasks[task_id].push_back({task_id, processor_id, start_time});
        result_processors[processor_id].push_back({task_id, start_time});
    }

    for (auto task_id : task_ids) {
        if (result_task_ids.find(task_id) == result_task_ids.end()) {
            std::cout << "Task " << task_id << " not assigned!!!" << std::endl;
            return 0;
        }
    }
    if (task_ids.size() < result_task_ids.size()) {
        for (auto task_id : result_task_ids) {
            if (task_ids.find(task_id) == task_ids.end()) {
                std::cout << "Invalid task id " << task_id << "!!!" << std::endl;
                return 0;
            }
        }
    }

    for (auto processor_id : result_processor_ids) {
        if (processor_ids.find(processor_id) == processor_ids.end()) {
            std::cout << "Invalid processor id " << processor_id << "!!!" << std::endl;
            return 0;
        }
    }

    std::unordered_map<int, double> task_timeout_response;
    std::unordered_map<int, double> task_efficiency;
    double sumMinTime = 0.0, sumProcessTime = 0.0;
    for (auto& p : result_tasks) {
        auto& assignments = p.second;
        if (assignments.size() > 5) {
            std::cout << "Task " << p.first << " assigned " << assignments.size() << " times!!!" << std::endl;
            return 0;
        }

        if (assignments.size() > 1) {
            std::sort(assignments.begin(), assignments.end(),
                      [](const Assignment& lhs, const Assignment& rhs)
                      { return lhs.start_time < rhs.start_time; });
        }

        int task_id = assignments[0].task_id;
        int processor_id = assignments[0].processor_id;
        int start_time = assignments[0].start_time;
        const Task& task = tasks[task_id];

        if (start_time < task.time) {
            std::cout << "Invalid earlier assignment time " << start_time << " for task " << task_id << "!!!" << std::endl;
            return 0;
        }
        if (start_time > task.time + task.sla_time) {
            //std::cout << "Assignment time " << start_time << " for task " << task_id << " beyond normal response range!" << std::endl;
            task_timeout_response[task_id] = (start_time - (task.time + task.sla_time)) / double(task.sla_time);
        }
        if (assignments.size() > 1) {
            for (std::size_t i = 1; i < assignments.size(); ++i) {
                int cur_processor_id = assignments[i].processor_id;
                int cur_start_time = assignments[i].start_time;
                if (processor_id == cur_processor_id) {
                    std::cout << "Task " << task_id << " continuously assigned to the same processor " << processor_id << "!!!" << std::endl;
                    return 0;
                }
                if (start_time == cur_start_time) {
                    std::cout << "Assigned task " << task_id << " should stay with processor " << processor_id << " at least 1 minute!!!" << std::endl;
                    return 0;
                }
                if (start_time + pt_matrix[processor_id][tasks[task_id].question_id] - 1 < cur_start_time) {
                    std::cout << "Task " << task_id << " finished before reassigned to processor " << cur_processor_id << "!!!" << std::endl;
                    return 0;
                }
                processor_id = cur_processor_id;
                start_time = cur_start_time;
            }
        }
        const auto& last_assignment = assignments.back();
        processor_id = last_assignment.processor_id;
        start_time = last_assignment.start_time;
        int process_time = pt_matrix[processor_id][tasks[task_id].question_id];
        int finish_time = start_time + process_time - 1;
        task_efficiency[task_id] = (double)process_time / (finish_time - task.time + 1);
        sumMinTime += minTime[task.question_id];
        sumProcessTime += (finish_time - task.time + 1);
    }
    double sum_timeout_response{0}, sum_efficiency{0};
    for (auto task_id : task_ids) {
        sum_timeout_response += task_timeout_response[task_id];
        sum_efficiency += task_efficiency[task_id];
    }
    double avg_timeout_response = sum_timeout_response / task_ids.size();
    double avg_efficiency = sum_efficiency / task_ids.size();

    std::unordered_map<int, double> processor_loads;
    for (const auto& p : result_processors) {
        std::unordered_map<int, int> time_cnt;
        int total_cnt{0};
        int processor_id = p.first;
        for (const auto& pp : p.second) {
            int task_id = pp.first;
            int question_id = tasks[task_id].question_id;
            int start_time = pp.second;
            if (result_tasks[task_id].size() > 1) {
                const auto& assignments = result_tasks[task_id];
                bool counted = false;
                for (std::size_t i = 0; i < assignments.size() - 1; ++i) {
                    if (start_time == assignments[i].start_time) {
                        for (int j = start_time; j < assignments[i+1].start_time; ++j) {
                            time_cnt[j] += 1;
                            total_cnt += 1;
                            if (time_cnt[j] > 3) {
                                std::cout << "More than 3 tasks for processor " << processor_id
                                          << " at time " << j << "!!!" << std::endl;
                                return 0;
                            }
                        }
                        counted = true;
                        break;
                    }
                }
                if (counted) continue;
            }
            for (int i = start_time; i < start_time + pt_matrix[processor_id][question_id]; ++i) {
                time_cnt[i] += 1;
                total_cnt += 1;
                if (time_cnt[i] > 3) {
                    std::cout << "More than 3 tasks for processor " << processor_id << " at time " << i << "!!!" << std::endl;
                    return 0;
                }
            }
        }
        processor_loads[processor_id] = total_cnt / (60.0 * 8 * 3);
    }

    double total_loads{0};
    for (auto processor_id : processor_ids) {
        printf("%d %f\n", processor_id, processor_loads[processor_id]);
        total_loads += processor_loads[processor_id];
    }
    double avg_loads = total_loads / processor_ids.size();
    double sum_loads{0};
    for (auto processor_id : processor_ids) {
        sum_loads += std::pow(processor_loads[processor_id] - avg_loads, 2);
    }
    double delta_loads = std::sqrt(sum_loads / processor_ids.size());

    int a{99}, b{99}, c{1000};
    double score = c * sumMinTime / sumProcessTime - (a * avg_timeout_response + b * delta_loads);

    printf("ave = %f\n", avg_loads);
    std::cout << "R: " << sumMinTime / sumProcessTime << "; M: " << avg_timeout_response << "; deltaL: " << delta_loads << std::endl;
    std::cout << "score: " << score << std::endl;

    return 0;
}
