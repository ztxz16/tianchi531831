//
// Created by huangyuyang on 10/21/20.
//

#include <map>
#include <vector>
#include <set>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <algorithm>

using namespace std;

string resultFile = "r.csv";
string taskFile = "/tcdata/sic_semi_work_order.csv";
string matrixFile = "/tcdata/sic_semi_process_time_matrix.csv";

void PrintResult();

clock_t startClock;
void CheckTime() {
    if ((double)(clock() - startClock) / CLOCKS_PER_SEC > 60 * 40) {
        PrintResult();
        exit(0);
    }
}

map <int, map <int, int> > matrix; // (pid, tid, time)
set <int> pset;

vector <int> split(string s, char c) {
    vector <int> ret;
    string cur = "";
    for (int i = 0; i < s.size(); i++) {
        if (s[i] == c) {
            ret.push_back(atoi(cur.c_str()));
            cur = "";
        } else {
            cur += s[i];
        }
    }
    if (cur != "") {
        ret.push_back(atoi(cur.c_str()));
    }
    return ret;
}
char buffer[1000005];
void ReadMatrix() {
    FILE *fi = fopen(matrixFile.c_str(), "r");
    int ret = fscanf(fi, "%s", buffer);
    vector <int> head = split(buffer, ',');
    while (fscanf(fi, "%s", buffer) == 1) {
        vector <int> p = split(buffer, ',');
        for (int i = 1; i < p.size(); i++) {
            matrix[p[0]][head[i]] = p[i];
        }
        pset.insert(p[0]);
    }
    fclose(fi);
}

struct PTask {
    int pid, tid, l, r;

    PTask (int pid, int tid, int l, int r) :
            pid(pid), tid(tid), l(l), r(r) {}
};

struct PTaskCompre {
    bool operator () (PTask *a, PTask *b) {
        vector <int> ca = vector <int> {a->r - a->l, a->l, a->tid, a->pid};
        vector <int> cb = vector <int> {b->r - b->l, b->l, b->tid, b->pid};
        return ca > cb;
    }
};

struct Task {
    int minTime;
    int finishTime, spendTime, acTime;
    int id, startTime, type, sla;
    set <PTask*, PTaskCompre> tasks;

    bool Read(FILE *fi) {
        char buffer[1005];
        if (fscanf(fi, "%s", buffer) != 1) {
            return false;
        }
        int pos = 0;
        while (buffer[pos] < '0' || buffer[pos] > '9') {
            pos++;
        }
        bool ok = sscanf(buffer + pos, "%d,%d,%d,%d", &id, &startTime, &type, &sla) == 4;
        if (ok) {
            minTime = 1e9;
            for (int p : pset) {
                minTime = min(minTime, matrix[p][type]);
            }
            return true;
        } else {
            return false;
        }
    }

    void Update() {
        finishTime = -1;
        acTime = 1e9;
        for (PTask *task : tasks) {
            finishTime = max(finishTime, task->r);
            acTime = min(acTime, task->l);
        }
        spendTime = finishTime - startTime;
    }

    void Insert(PTask *task) {
        tasks.insert(task);
        Update();
    }

    void Remove(PTask *task) {
        tasks.erase(task);
        Update();
    }
};
map <int, Task> tasks; // (tid -> taskList)
void ReadTask() {
    FILE *fi = fopen(taskFile.c_str(), "r");
    Task cur;
    while (cur.Read(fi)) {
        tasks[cur.id] = cur;
    }
    fclose(fi);
}

struct People {
    int times; // cost times
    set <PTask*, PTaskCompre> tasks;
    map <int, int> payloads;

    bool CanInsert(PTask *task) {
        int l = task->l;
        int r = task->r;
        for (int i = l; i < r; i++) {
            if (payloads[i] > 2) {
                return false;
            }
        }
        return true;
    }

    void Add(int l, int r, int v) {
        for (int i = l; i < r; i++) {
            payloads[i] += v;
        }
    }

    void Insert(PTask *task) {
        int l = task->l;
        int r = task->r;
        for (int i = l; i < r; i++) {
            payloads[i]++;
        }
        tasks.insert(task);
        times += (r - l);
    }

    void Remove(PTask *task) {
        int l = task->l;
        int r = task->r;
        for (int i = l; i < r; i++) {
            payloads[i]--;
        }
        tasks.erase(task);
        times -= (r - l);
    }
};
map <int, People> people;

bool IsFirstTask(PTask *task) {
    return tasks[task->tid].acTime == task->l;
}

bool IsFinishTask(PTask *task) {
    return tasks[task->tid].finishTime == task->r;
}

bool CanInsertPTask(PTask *task) {
    return people[task->pid].CanInsert(task);
}

void InsertPTask(PTask *task) {
    people[task->pid].Insert(task);
    tasks[task->tid].Insert(task);
}

void RemovePTask(PTask *task) {
    people[task->pid].Remove(task);
    tasks[task->tid].Remove(task);
}

void ReadResult() {
    FILE *fi = fopen(resultFile.c_str(), "r");
    map <int, vector <pair <int, int> > > taskList; // taskId -> {time, pid}
    int tid, pid, time;
    while (fscanf(fi, "%d,%d,%d", &tid, &pid, &time) == 3) {
        taskList[tid].push_back(make_pair(time, pid));
    }
    for (auto &it : taskList) {
        sort(it.second.begin(), it.second.end());
        for (int i = 0; i < it.second.size(); i++) {
            int pid = it.second[i].second, l = it.second[i].first;
            int r = i < it.second.size() - 1 ? it.second[i + 1].first : l + matrix[pid][tasks[it.first].type];
            PTask *task = new PTask(pid, it.first, l, r);
            InsertPTask(task);
        }
    }
    fclose(fi);
}

struct Score {
    double R, M, L, score;

    void print() {
        printf("R = %f, M = %f, L = %f, score = %f\n", R, M, L, score);
    }
};

Score Calc() {
    // R, M
    double sumMin = 0;
    double sumSpend = 0;
    double sumDelay = 0;
    for (auto &it : tasks) {
        sumMin += it.second.minTime;
        sumSpend += it.second.spendTime;

        if (it.second.acTime > it.second.startTime + it.second.sla) {
            sumDelay += (double)(it.second.acTime - it.second.startTime - it.second.sla) / it.second.sla;
        }
    }
    double R = sumMin / sumSpend;
    double M = sumDelay / tasks.size();

    // L
    double sum = 0.0;
    for (int pid : pset) {
        double cur = (double)people[pid].times / (60.0 * 8 * 3);
        sum += cur;
    }
    double ave = sum / pset.size();
    sum = 0;
    for (int pid : pset) {
        double cur = (double)people[pid].times / (60.0 * 8 * 3);
        sum += (ave - cur) * (ave - cur);
    }
    double L = sqrt(sum / pset.size());
    double score = R * 1000 - 99 * (M + L);

    Score ret;
    ret.R = R;
    ret.M = M;
    ret.L = L;
    ret.score = score;
    return ret;
}

void PrintResult() {
    FILE *fo = fopen("r.csv", "w");
    for (int i : pset) {
        set <PTask*, PTaskCompre> tasks = people[i].tasks;
        for (PTask *task : tasks) {
            fprintf(fo, "%d,%d,%d\n", task->tid, i, task->l);
        }
    }
    fclose(fo);
}

bool OptSwap(double &curResult) {
    bool ret = false;
    for (int i : pset) {
        set <PTask*, PTaskCompre> tasks = people[i].tasks;
        for (auto &task : tasks) {
            // TODO: 任意情况都做调整
            if (!IsFinishTask(task)) {
                continue;
            }

            if (::tasks[task->tid].tasks.size() == 1) {
                int type = ::tasks[task->tid].type;
                int sel = -1, start = -1;
                for (int j : pset) {
                    if (i == j) {
                        continue;
                    }
                    int spend = matrix[j][type];
                    PTask *newTask = new PTask(j, task->tid, ::tasks[task->tid].startTime,
                                               ::tasks[task->tid].startTime + spend);
                    while (newTask->r < task->r) {
                        if (CanInsertPTask(newTask)) {
                            break;
                        }
                        newTask->l++;
                        newTask->r++;
                    }
                    if (newTask->r < task->r - 1) {
                        InsertPTask(newTask);
                        double cur = Calc().score;
                        if (cur > curResult) {
                            curResult = cur;
                            sel = j;
                            start = newTask->l;
                        }
                        RemovePTask(newTask);
                    }
                    delete newTask;
                }
                if (sel != -1) {
                    ret = true;
                    PTask *newTask = new PTask(sel, task->tid, start, start + matrix[sel][type]);
                    RemovePTask(task);
                    InsertPTask(newTask);
                }
            } else {
                /*int type = ::tasks[task->tid].type;
                int sel = -1, start = -1;
                for (int j : pset) {
                    if (i == j) {
                        continue;
                    }
                    int spend = matrix[j][type];
                    PTask *newTask = new PTask(j, task->tid, ::tasks[task->tid].startTime,
                                               ::tasks[task->tid].startTime + spend);
                    while (newTask->r < task->r) {
                        if (CanInsertPTask(newTask)) {
                            break;
                        }
                        newTask->l++;
                        newTask->r++;
                    }
                    if (newTask->r < task->r - 1 && newTask->l <= task->l) {
                        InsertPTask(newTask);
                        double cur = Calc().score;
                        if (cur > curResult) {
                            curResult = cur;
                            sel = j;
                            start = newTask->l;
                        }
                        RemovePTask(newTask);
                    }
                    delete newTask;
                }
                if (sel != -1) {
                    PTask *newTask = new PTask(sel, task->tid, start, start + matrix[sel][type]);
                    set <PTask*, PTaskCompre> oriTasks = ::tasks[task->tid].tasks;
                    for (auto task : oriTasks) {
                        RemovePTask(task);
                    }
                    InsertPTask(newTask);
                    for (auto task : oriTasks) {
                        if (task->l >= newTask->l) {
                            delete task;
                        } else if (task->r > newTask->l) {
                            PTask *curNewTask = new PTask(task->pid, task->tid, task->l, newTask->l);
                            InsertPTask(curNewTask);
                            delete task;
                        } else {
                            InsertPTask(task);
                        }
                    }
                }*/
            }
        }
    }
    return ret;
}

bool OptLeftMove(double &curResult) {
    bool ret = false;
    for (int i : pset) {
        set <PTask*, PTaskCompre> tasks = people[i].tasks;
        for (auto &task : tasks) {
            if (!IsFinishTask(task)) {
                continue;
            }
            RemovePTask(task);
            int spend = task->r - task->l;
            for (int l = ::tasks[task->tid].startTime; l <= task->l; l++) {
                int r = l + spend;
                PTask *newTask = new PTask(i, task->tid, l, r);
                if (CanInsertPTask(newTask)) {
                    if (l < task->l) {
                        ret = true;
                    }
                    set <PTask*, PTaskCompre> oriTasks = ::tasks[task->tid].tasks;
                    for (auto ori : oriTasks) {
                        if (ori->l >= newTask->l) {
                            RemovePTask(ori);
                            delete ori;
                        } else if (ori->r > newTask->l) {
                            RemovePTask(ori);
                            ori->r = newTask->l;
                            InsertPTask(ori);
                        }
                    }

                    InsertPTask(newTask);
                    delete task;
                    break;
                } else {
                    delete newTask;
                }
            }
        }
    }
    return ret;
}

void OptSimpleSwap(double &curResult) {
    curResult = Calc().score;
    for (int i : pset) {
        set <PTask*, PTaskCompre> tasks = people[i].tasks;
        for (auto &task : tasks) {
            set <int> curP;
            for (auto &it : ::tasks[task->tid].tasks) {
                curP.insert(it->pid);
            }
            int type = ::tasks[task->tid].type;
            int sel = -1;
            bool isFinish = IsFinishTask(task);
            RemovePTask(task);
            for (int j : pset) {
                // TODO: 合并同一个PID的连续任务
                if (curP.find(j) != curP.end()) {
                    continue;
                }
                PTask *newTask = new PTask(j, task->tid, task->l, task->r);
                if (isFinish) {
                    if (matrix[i][type] + 5 < matrix[j][type]) {
                        continue;
                    }
                    newTask->r = newTask->l + matrix[j][type];
                }
                if (CanInsertPTask(newTask)) {
                    InsertPTask(newTask);
                    Score score = Calc();
                    if (score.score > curResult) {
                        curResult = score.score;
                        sel = j;
                    }
                    RemovePTask(newTask);
                }
                delete newTask;
            }
            if (sel != -1) {
                PTask *newTask = new PTask(sel, task->tid, task->l, task->r);
                if (isFinish) {
                    newTask->r = newTask->l + matrix[sel][type];
                }
                InsertPTask(newTask);
                delete task;
            } else {
                InsertPTask(task);
            }
        }
    }
}

void OptSimpleSwapAndSplit(double &curResult) {
    curResult = Calc().score;
    for (int i : pset) {
        CheckTime();
        set <PTask*, PTaskCompre> tasks = people[i].tasks;
        for (auto &task : tasks) {
            set <int> curP;
            for (auto &it : ::tasks[task->tid].tasks) {
                curP.insert(it->pid);
            }
            int type = ::tasks[task->tid].type;
            int sel = -1, start = -1, end = -1;
            if (IsFinishTask(task) || ::tasks[task->tid].tasks.size() >= 5) {
                continue;
            }
            RemovePTask(task);
            for (int j : pset) {
                if (curP.find(j) != curP.end()) {
                    continue;
                }
                PTask *newTask = new PTask(j, task->tid, task->l, task->r);
                while (newTask->r > newTask->l) {
                    if (CanInsertPTask(newTask)) {
                        PTask *oldTask = new PTask(i, task->tid, newTask->r, task->r);
                        InsertPTask(newTask);
                        InsertPTask(oldTask);
                        Score score = Calc();
                        if (score.score > curResult) {
                            curResult = score.score;
                            sel = j;
                            start = newTask->l;
                            end = newTask->r;

                            printf("%f\n", curResult);
                        }
                        RemovePTask(oldTask);
                        RemovePTask(newTask);
                        delete oldTask;

                        //break;
                    }
                    newTask->r--;
                }

                newTask->l = task->l;
                newTask->r = task->r;
                while (newTask->r > newTask->l) {
                    if (CanInsertPTask(newTask)) {
                        PTask *oldTask = new PTask(i, task->tid, task->l, newTask->l);
                        InsertPTask(newTask);
                        InsertPTask(oldTask);
                        Score score = Calc();
                        if (score.score > curResult) {
                            curResult = score.score;
                            sel = j;
                            start = newTask->l;
                            end = newTask->r;

                            printf("%f\n", curResult);
                        }
                        RemovePTask(oldTask);
                        RemovePTask(newTask);
                        delete oldTask;

                        //break;
                    }
                    newTask->l++;
                }

                delete newTask;
            }
            if (sel != -1) {
                PTask *newTask = new PTask(sel, task->tid, start, end);
                InsertPTask(newTask);
                if (end < task->r) {
                    PTask *oldTask = new PTask(i, task->tid, end, task->r);
                    InsertPTask(oldTask);
                } else if (start > task->l) {
                    PTask *oldTask = new PTask(i, task->tid, task->l, start);
                    InsertPTask(oldTask);
                }
                delete task;
            } else {
                InsertPTask(task);
            }
        }
    }
}

bool OptSLA(double &curResult) {
    vector <pair <pair <int, int>, int> > taskList;
    for (auto &it : tasks) {
        taskList.push_back(make_pair(make_pair(it.second.acTime <= it.second.startTime + it.second.sla, it.second.sla), it.first));
    }
    sort(taskList.begin(), taskList.end());
    bool ret = false;
    for (int i = 0; i < taskList.size(); i++) {
        auto it = *tasks.find(taskList[i].second);
        auto task = it.second;
        if (task.tasks.size() >= 5) {
            continue;
        }

        set <int> curP;
        for (auto p : it.second.tasks) {
            curP.insert(p->pid);
        }
        int type = task.acTime;
        int sel = -1;
        int minL = task.acTime;
        int start = -1;
        for (int j : pset) {
            if (curP.find(j) != curP.end()) {
                continue;
            }
            PTask *newTask = new PTask(j, it.first, task.startTime, task.acTime);
            while (newTask->l < minL) {
                if (CanInsertPTask(newTask)) {
                    break;
                }
                newTask->l++;
            }
            /*
            if (newTask->l < minL) {
                sel = j;
                minL = newTask->l;
            }
            delete newTask;
             */
            if (newTask->l < minL) {
                InsertPTask(newTask);
                Score score = Calc();
                if (score.score > curResult) {
                    curResult = score.score;
                    sel = j;
                    start = newTask->l;
                }
                RemovePTask(newTask);
            }
        }
        if (sel != -1) {
            PTask *newTask = new PTask(sel, it.first, start, task.acTime);
            InsertPTask(newTask);
            ret = true;
            i--;
        }
    }
    return ret;
}

int main() {
    startClock = clock();

    ReadMatrix();
    ReadTask();
    ReadResult();

    Score s = Calc();
    s.print();

    double curResult = s.score;
/*
    while (true) {
        double lastR = Calc().R;
        OptLeftMove(curResult);
        Calc().print();
        OptSwap(curResult);
        Calc().print();
        PrintResult();
        if (Calc().R < lastR + 1e-9) {
            break;
        }
    }
    return 0;
*/
    while (true) {
        double lastResult = curResult;
        while (true) {
            double lastResult = curResult;
            //OptSimpleSwap(curResult);
            OptLeftMove(curResult);
            OptSwap(curResult);

            if (curResult < lastResult + 1e-9) {
                break;
            }
        }

        printf("opt 0: ");
        Calc().print();
        while (OptSLA(curResult));
        printf("opt 1: ");
        Calc().print();

        while (true) {
            double lastResult = curResult;
            OptSimpleSwap(curResult);
            if (curResult < lastResult + 1e-9) {
                break;
            }
        }

        PrintResult();
        while (true) {
            double lastResult = curResult;
            OptSimpleSwapAndSplit(curResult);
            if (curResult < lastResult + 1e-9) {
                break;
            }
        }

        PrintResult();
        if (curResult < lastResult + 1e-9) {
            break;
        }

        printf("opt 2: ");
        Calc().print();
    }

    PrintResult();
    return 0;
}