#include <set>
#include <map>
#include <vector>
#include <string>
#include <algorithm>

using namespace std;

const int MCFLIM = 1000005;
struct MCF {
    int e[MCFLIM], succ[MCFLIM], last[MCFLIM], val[MCFLIM], cost[MCFLIM], sum;
    int dis[MCFLIM], pre[MCFLIM], queue[MCFLIM];

    void init(int n) {
        sum = 1;
        for (int i = 0; i <= n; i++) {
            last[i] = 0;
        }
    }

    void add(int a, int b, int c, int d) {
        e[++sum] = b, succ[sum] = last[a], last[a] = sum;
        e[++sum] = a, succ[sum] = last[b], last[b] = sum;
        val[sum - 1] = c, val[sum] = 0;
        cost[sum - 1] = d, cost[sum] = -d;
    }

    void spfa(int s, int t) {
        int h, tt, now, x;
        for (int i = 1; i <= t; i++) {
            dis[i] = 1e9;
        }
        dis[s] = 0;
        queue[1] = s;
        for (h = tt = 1; h <= tt; h++) {
            now = queue[h];
            for (x = last[now]; x != 0; x = succ[x])
                if (val[x] && dis[now] + cost[x] < dis[e[x]])
                    dis[e[x]] = dis[now] + cost[x], queue[++tt] = e[x], pre[e[x]] = x;
        }
    }

    void work(int s, int t, int &ans1, int &ans2) {
        int now, x = 1e9;
        for (now = t; now != s; now = e[pre[now] ^ 1])
            x = min(x, val[pre[now]]);
        ans1 += x;
        for (now = t; now != s; now = e[pre[now] ^ 1])
            val[pre[now]] -= x, val[pre[now] ^ 1] += x, ans2 += x * cost[pre[now]];
    }

    int getAns(int s, int t) {
        int ans1 = 0, ans2 = 0;
        while (1) {
            spfa(s, t);
            if (dis[t] == (int)1e9)
                break;
            work(s, t, ans1, ans2);
        }
        return ans2;
    }

    void getPairs(int s, int mid, vector <pair <int, int> > &pairs) {
        pairs.clear();
        for (int i = 1; i <= mid; i++) {
            for (int x = last[i]; x != 0; x = succ[x]) {
                if (e[x] != s && val[x] == 0 && cost[x] < 9999) {
                    pairs.push_back(make_pair(i, e[x] - mid));
                }
            }
        }
    }
};

string taskFile = "/tcdata/sic_semi_work_order.csv";
string matrixFile = "/tcdata/sic_semi_process_time_matrix.csv";

map <pair <int, int>, int> matrix; // (pid, tid, time)
vector <pair <int, int> > plist; // (pid, finishTime)

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
            matrix[make_pair(p[0], head[i])] = p[i];
        }
        plist.push_back(make_pair(p[0], 0));
        plist.push_back(make_pair(p[0], 0));
        plist.push_back(make_pair(p[0], 0));
    }
    fclose(fi);
}

struct Task {
    int id, startTime, type, sla;
    int disTimes = 0, belong = -1, finishTime = -1;

    bool Read(FILE *fi) {
        char buffer[1005];
        if (fscanf(fi, "%s", buffer) != 1) {
            return false;
        }
        int pos = 0;
        while (buffer[pos] < '0' || buffer[pos] > '9') {
            pos++;
        }
        return sscanf(buffer + pos, "%d,%d,%d,%d", &id, &startTime, &type, &sla) == 4;
    }

    bool operator < (const Task &b) {
        return startTime + sla < b.startTime + b.sla;
    }
};
map <int, vector <Task> > tasks; // (time -> taskList)

struct TaskCompre {
    bool operator () (Task *a, Task *b) {
        return a->id > b->id;
        vector <int> va = vector <int> {a->startTime, a->id};
        vector <int> vb = vector <int> {b->startTime, b->id};
        return va > vb;
    }
};

void ReadTask() {
    FILE *fi = fopen(taskFile.c_str(), "r");
    Task cur;
    while (cur.Read(fi)) {
        tasks[cur.startTime + cur.sla].push_back(cur);
    }
    fclose(fi);
}

MCF mcf;

int main() {
    ReadMatrix();
    ReadTask();
    FILE *fo = fopen("result.csv", "w");
    set <Task*, TaskCompre> curTasks;
    double rsum = 0, rcnt = 0;
    int lastTime = 0;
    for (int tt = 0; tt <= 1440 || curTasks.size() > 0; tt++) {
        vector <int> curP;
        for (int i = 0; i < plist.size(); i++) {
            if (plist[i].second <= tt || true) {
                curP.push_back(i);
            }
        }
        if (tasks.find(tt) != tasks.end()) {
            for (auto &task : tasks[tt]) {
                curTasks.insert(&task);
            }
        }

        //printf("%d %d %d\n", tt, (int)curTasks.size(), (int)curP.size());
        if ((curTasks.size() >= 280 && curP.size() >= 270) ||
            tt > 1440) {
            vector <Task*> vtasks;
            for (auto task : curTasks) {
                vtasks.push_back(task);
            }
            int n = curP.size(), m = vtasks.size(), s = n + m + 1, t = s + 1;
            if (n == 0 || m == 0) {
                continue;
            }
            vector <pair <int, int> > pairs;
            mcf.init(t);
            for (int a = 1; a <= n; a++) {
                mcf.add(s, a, 1, 0);
            }
            for (int a = 1; a <= m; a++) {
                mcf.add(n + a, t, 1, 0);
            }
            for (int a = 1; a <= n; a++) {
                int pid = plist[curP[a - 1]].first;
                for (int b = 1; b <= m; b++) {
                    int startTime = max(plist[curP[a - 1]].second, vtasks[b - 1]->startTime);
                    if (matrix[make_pair(pid, vtasks[b - 1]->type)] < 9999) {
	                    mcf.add(a, n + b, 1, startTime + matrix[make_pair(pid, vtasks[b - 1]->type)]);
                    }
                }
            }
            mcf.getAns(s, t);
            mcf.getPairs(s, n, pairs);

            set <int> emptyP;
            for (int p : curP) {
                emptyP.insert(p);
            }
            for (auto it : pairs) {
                emptyP.erase(curP[it.first - 1]);
            }
            int nextFinishTime = 1e9;
            for (int i = 0; i < plist.size(); i++) {
                if (plist[i].second > tt) {
                    nextFinishTime = min(nextFinishTime, plist[i].second);
                }
            }
            for (auto it : pairs) {
                Task *task = vtasks[it.second - 1];
                int pindex = curP[it.first - 1];
                int pid = plist[pindex].first;
                int startTime = max(plist[pindex].second, task->startTime);
                nextFinishTime = min(nextFinishTime, startTime + matrix[make_pair(pid, task->type)]);
            }

            for (auto it : pairs) {
                Task *task = vtasks[it.second - 1];
                int pindex = curP[it.first - 1];
                int pid = plist[pindex].first;
                int startTime = max(plist[pindex].second, task->startTime);
                if (startTime >= nextFinishTime) {
                    continue;
                }
                double curR = -1;

                if (emptyP.size() > 0 && startTime > task->startTime + task->sla) {
                    int slaTime = task->startTime + task->sla;
                    int sel = -1, time = startTime;
                    for (auto p : emptyP) {
                        if ((time > slaTime && plist[p].second < time) ||
                            (time <= slaTime && plist[p].second <= slaTime && plist[p].second > time)) {
                            sel = p;
                            time = plist[p].second;
                        }
                    }
                    if (sel != -1) {
                        curR = (double)max(0, max(task->startTime, plist[sel].second) - task->startTime - task->sla) / task->sla;
                        //fprintf(fo, "%d,%d,%d\n", task->id, plist[sel].first, max(task->startTime, plist[sel].second));

                        //printf("pid: %d, gap: (%d -> %d -> %d)\n", plist[sel].first, plist[sel].second, max(task->startTime, plist[sel].second), startTime);
                        plist[sel].second = startTime;
                        emptyP.erase(sel);
                    }
                }

                fprintf(fo, "%d,%d,%d\n", task->id, pid, startTime);
                if (curR < 0) {
                    curR = (double)max(0, startTime - task->startTime - task->sla) / task->sla;
                }
                rsum += curR;
                rcnt++;
                //printf("rsum = %f, rcnt = %f, rave = %f\n", rsum, rcnt, rsum / rcnt);
                //printf("pid: %d, gap: (%d -> %d -> %d)\n", plist[pindex].first, plist[pindex].second, startTime, startTime + matrix[make_pair(pid, task->type)]);
                plist[pindex].second = startTime + matrix[make_pair(pid, task->type)];
                curTasks.erase(task);
            }
        }
    }

    printf("rsum = %f, rcnt = %f, rave = %f\n", rsum, rcnt, rsum / rcnt);
    fclose(fo);
    return 0;
}