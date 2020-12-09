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
#include <ctime>

using namespace std;

namespace Worker {
    string resultFile = "r.csv";
    string taskFile = "/tcdata/sic_semi_work_order.csv";
    string matrixFile = "/tcdata/sic_semi_process_time_matrix.csv";

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

    void PrintResult(string fileName) {
        FILE *fo = fopen(fileName.c_str(), "w");
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

                if (Worker::tasks[task->tid].tasks.size() == 1) {
                    int type = Worker::tasks[task->tid].type;
                    int sel = -1, start = -1;
                    for (int j : pset) {
                        if (i == j) {
                            continue;
                        }
                        int spend = matrix[j][type];
                        PTask *newTask = new PTask(j, task->tid, Worker::tasks[task->tid].startTime,
                                                   Worker::tasks[task->tid].startTime + spend);
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
                for (int l = Worker::tasks[task->tid].startTime; l <= task->l; l++) {
                    int r = l + spend;
                    PTask *newTask = new PTask(i, task->tid, l, r);
                    if (CanInsertPTask(newTask)) {
                        if (l < task->l) {
                            ret = true;
                        }
                        set <PTask*, PTaskCompre> oriTasks = Worker::tasks[task->tid].tasks;
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

    double redo() {
        matrix.clear();
        pset.clear();
        people.clear();
        tasks.clear();
        ReadMatrix();
        ReadTask();
        ReadResult();

        while (true) {
            double curResult;
            double lastR = Calc().R;
            OptLeftMove(curResult);
            Calc().print();
            OptSwap(curResult);
            Calc().print();
            PrintResult("r.csv");
            if (Calc().R < lastR + 1e-9) {
                break;
            }
        }

        return Calc().R;
    }
}

string resultFile = "result.csv";
string taskFile = "/tcdata/sic_semi_work_order.csv";
string matrixFile = "/tcdata/sic_semi_process_time_matrix.csv";

struct Task;
map <int, Task> tasks;
set <int> pset;
map <int, map <int, int> > matrix; // (pid, tid, time)

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

struct TaskCompare {
    bool operator () (Task *a, Task *b);
};
struct People {
    int pid;
    set <Task*, TaskCompare> tasks;
    int finishTime;

    People (int pid) : pid(pid), finishTime(0) {}

    double GetSum();
};
vector <People> peoples;
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
        peoples.push_back(People(p[0]));
        peoples.push_back(People(p[0]));
        peoples.push_back(People(p[0]));
    }
    fclose(fi);
}

struct Task {
    int minTime;
    int id, startTime, type, sla;

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
};

bool TaskCompare::operator()(Task *a, Task *b) {
    return vector<int> {a->startTime, a->minTime, a->id} < vector<int> {b->startTime, b->minTime, b->id};
}

double startSum = 0.0f, minSum = 0.0f;
void ReadTask() {
    FILE *fi = fopen(taskFile.c_str(), "r");
    Task cur;
    while (cur.Read(fi)) {
        tasks[cur.id] = cur;
        startSum += cur.startTime;
        minSum += cur.minTime;
    }
    fclose(fi);
}

double People::GetSum()  {
    int cur = 0;
    double sum = 0;
    int cnt = 0;
    /*for (auto task : tasks) {
        int startTime = max(cur, task->startTime);
        int endTime = startTime + matrix[pid][task->type];
        cur = endTime;
        sum += cur;
    }*/
    multiset <int> s;
    auto tt = tasks;
    while (cnt < tasks.size()) {
        while (tt.size() > 0 && (*tt.begin())->startTime <= cur) {
            s.insert(matrix[pid][(*tt.begin())->type]);
            tt.erase(tt.begin());
        }
        if (s.size() > 0) {
            cnt++;
            cur += *s.begin();
            s.erase(s.begin());
            sum += cur;
        } else {
            cur = (*tt.begin())->startTime;
        }
    }
    return sum;
}

bool PTaskCMP(PTask *a, PTask *b) {
    return vector <int> {a->l, a->r, a->pid, a->tid} <
           vector <int> {b->l, b->r, b->pid, b->tid};
}

void ReadResult(string resultFile) {
    FILE *fi = fopen(resultFile.c_str(), "r");
    vector <PTask*> tasks;
    int tid, pid, time;
    while (fscanf(fi, "%d,%d,%d", &tid, &pid, &time) == 3) {
        tasks.push_back(new PTask(pid, tid, time, time + matrix[pid][::tasks[tid].type]));
    }
    sort(tasks.begin(), tasks.end(), PTaskCMP);
    for (auto &task : tasks) {
        int pid = task->pid;
        int sel = -1, maxTime = -1;
        for (int j = 0; j < peoples.size(); j++) {
            if (peoples[j].pid == pid && peoples[j].finishTime <= task->l && peoples[j].finishTime > maxTime) {
                maxTime = peoples[j].finishTime;
                sel = j;
            }
        }
        peoples[sel].tasks.insert(&::tasks[task->tid]);
        peoples[sel].finishTime = task->r;
    }

    fclose(fi);
}

double Calc() {
    double finishSum = 0.0;
    for (int i = 0; i < peoples.size(); i++) {
        finishSum += peoples[i].GetSum();
    }
    return minSum / (finishSum - startSum);
}

void PrintResult(string fileName) {
    FILE *fo = fopen(fileName.c_str(), "w");
    for (int i = 0; i < peoples.size(); i++) {
        multiset <pair <int, int> > s;
        auto tt = peoples[i].tasks;
        int cnt = 0;
        int cur = 0;
        int pid = peoples[i].pid;
        while (cnt < peoples[i].tasks.size()) {
            while (tt.size() > 0 && (*tt.begin())->startTime <= cur) {
                s.insert(make_pair(matrix[pid][(*tt.begin())->type], (*tt.begin())->id));
                tt.erase(tt.begin());
            }
            if (s.size() > 0) {
                fprintf(fo, "%d,%d,%d\n", s.begin()->second, pid, cur);
                cnt++;
                cur += s.begin()->first;
                s.erase(s.begin());
            } else {
                cur = (*tt.begin())->startTime;
            }
        }
        /*
        int cur = 0;
        for (auto task : peoples[i].tasks) {
            int startTime = max(cur, task->startTime);
            int endTime = startTime + matrix[peoples[i].pid][task->type];
            cur = endTime;
            fprintf(fo, "%d,%d,%d\n", task->id, peoples[i].pid, startTime);
        }
         */
    }
    fclose(fo);
}

clock_t startClock;
void CheckTime() {
    if ((double)(clock() - startClock) / CLOCKS_PER_SEC > 60 * 60 * 2 + 60 * 5) {
        exit(0);
    }
}

int main() {
    startClock = clock();
    ReadMatrix();
    ReadTask();
    ReadResult(resultFile);

    printf("R: %f\n", Calc());

    while (true) {
        int lastScore = Calc();
        if (true) {
            double bestScore = Calc();
            while (true) {
                double lastScore = bestScore;
                for (int i = 0; i < peoples.size(); i++) {
                    printf("i = %d\n", i);
                    CheckTime();
                    set<Task *, TaskCompare> cur = peoples[i].tasks;
                    for (Task *task : cur) {
                        double sumNow = minSum / Calc() + startSum - peoples[i].GetSum();
                        int sel = -1;
                        peoples[i].tasks.erase(task);
                        sumNow += peoples[i].GetSum();
                        for (int j = 0; j < peoples.size(); j++) {
                            if (matrix[peoples[j].pid][task->type] > 9999 || i == j) {
                                continue;
                            }
                            sumNow -= peoples[j].GetSum();
                            peoples[j].tasks.insert(task);
                            sumNow += peoples[j].GetSum();
                            double curScore = minSum / (sumNow - startSum);
                            if (curScore > bestScore) {
                                bestScore = curScore;
                                sel = j;

                                printf("opt 1 %f\n", bestScore);
                            }
                            sumNow -= peoples[j].GetSum();
                            peoples[j].tasks.erase(task);
                            sumNow += peoples[j].GetSum();

                            if (sel != -1) {
                                break;
                            }
                        }
                        if (sel != -1) {
                            peoples[sel].tasks.insert(task);
                        } else {
                            peoples[i].tasks.insert(task);
                        }
                    }
                }

                if (lastScore > bestScore - 1e-9) {
                    break;
                }
                printf("opt 1 %f\n", bestScore);
                PrintResult("r.csv");

                break;
            }
        }

        PrintResult("r.csv");
        Worker::redo();
        for (auto &it : peoples) {
            it.tasks.clear();
            it.finishTime = 0;
        }
        ReadResult("r.csv");

        if (true) {
            double bestScore = Calc();
            while (true) {
                double lastScore = bestScore;
                for (int i = 0; i < peoples.size(); i++) {
                    printf("i = %d\n", i);
                    CheckTime();
                    set<Task *, TaskCompare> cur = peoples[i].tasks;
                    for (Task *task : cur) {
                        double sumNow = minSum / Calc() + startSum - peoples[i].GetSum();
                        int sel = -1;
                        Task *another = nullptr;
                        peoples[i].tasks.erase(task);
                        sumNow += peoples[i].GetSum();
                        for (int j = 0; j < peoples.size(); j++) {
                            if (matrix[peoples[j].pid][task->type] > 9999 || i == j) {
                                continue;
                            }
                            sumNow -= peoples[j].GetSum();
                            peoples[j].tasks.insert(task);
                            sumNow += peoples[j].GetSum();
                            set<Task *, TaskCompare> jtasks = peoples[j].tasks;
                            for (Task *tt : jtasks) {
                                if (matrix[peoples[i].pid][tt->type] > 9999) {
                                    continue;
                                }

                                sumNow -= peoples[j].GetSum();
                                peoples[j].tasks.erase(tt);
                                sumNow += peoples[j].GetSum();
                                sumNow -= peoples[i].GetSum();
                                peoples[i].tasks.insert(tt);
                                sumNow += peoples[i].GetSum();

                                double curScore = minSum / (sumNow - startSum);
                                if (curScore > bestScore) {
                                    bestScore = curScore;
                                    sel = j;
                                    another = tt;

                                    printf("opt 2: %f\n", bestScore);
                                }

                                sumNow -= peoples[j].GetSum();
                                peoples[j].tasks.insert(tt);
                                sumNow += peoples[j].GetSum();
                                sumNow -= peoples[i].GetSum();
                                peoples[i].tasks.erase(tt);
                                sumNow += peoples[i].GetSum();
                            }

                            sumNow -= peoples[j].GetSum();
                            peoples[j].tasks.erase(task);
                            sumNow += peoples[j].GetSum();

                            if (sel != -1) {
                                break;
                            }
                        }

                        if (sel != -1) {
                            peoples[sel].tasks.insert(task);
                            peoples[sel].tasks.erase(another);
                            peoples[i].tasks.insert(another);
                        } else {
                            peoples[i].tasks.insert(task);
                        }
                    }
                }

                if (lastScore > bestScore - 1e-9) {
                    break;
                }
                printf("opt 2: %f\n", bestScore);
                PrintResult("r.csv");

                break;
            }
        }

        if (lastScore > Worker::redo() - 1e-9) {
            break;
        } else {
            for (auto &it : peoples) {
                it.tasks.clear();
                it.finishTime = 0;
            }
            ReadResult("r.csv");
        }
    }

    PrintResult("r.csv");
    Worker::redo();
    return 0;
}