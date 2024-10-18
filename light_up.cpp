#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <sstream>
#include <string>
#include <regex>
#include <utility>

#define MAXN 1024
#define MAXEDGES 10000
#define MAXNB 100

using namespace std;

int len = 1000, init_len = 1000;
const int repeat_times_bound = 1000;
int repeat_times = 0;
double sum = 0, sum2 = 0;
const double T_0 = 100, cooling_rate = 1 - 1e-4, T_min = 0.01;
double loss = 100000, loss_num = 100000;
int num = 0, sum_num = 0;
const double p_max = 6, p_min = 0, p_step = 1;
double p = 0, temp_p = 0, init_p = 0;
bool is_init = false; // change this to initialize the p
bool is_stable = false;

struct vertex
{
    int state;
    int neighbor[MAXNB];
    int size;
    vertex() : state(rand() % 2), size(0) {}
};

struct graph
{
    vertex vertices[MAXN];
    int edges[MAXEDGES][2];
    int vertices_size, edge_size;
    bool edge_matrix[MAXN][MAXN];
    double loss_1, loss_21, loss_31;

    void add_edge(int u, int v)
    {
        vertices[u].neighbor[vertices[u].size++] = v;
        vertices[v].neighbor[vertices[v].size++] = u;
        edges[edge_size][0] = u;
        edges[edge_size][1] = v;
        edge_size++;
        edge_matrix[u][v] = true;
        edge_matrix[v][u] = true;
    }

    graph(int n) : vertices_size(n), edge_size(0), loss_1(0), loss_21(0), loss_31(0)
    {
        fill(vertices, vertices + n, vertex());
        for (int i = 0; i < n; ++i)
            fill(edge_matrix[i], edge_matrix[i] + n, false);
        for (int i = 0; i < n; ++i)
        {
            for (int j = i + 1; j < n; ++j)
            {
                if (((double)rand() / RAND_MAX) < 0.01)
                {
                    add_edge(i, j);
                }
            }
        }
        for (int i = 0; i < vertices_size; ++i)
        {
            loss_1 += vertices[i].state;
        }
        for (int i = 0; i < edge_size; ++i)
        {
            loss_21 += double(vertices[edges[i][0]].state ^ vertices[edges[i][1]].state);
            loss_31 += double(vertices[edges[i][0]].state ^ vertices[edges[i][1]].state) * (1.0 / vertices[edges[i][0]].size + 1.0 / vertices[edges[i][1]].size);
        }
    }

    void light_up(int u)
    {
        vertices[u].state = 1 - vertices[u].state;
        loss_1 += 2 * vertices[u].state - 1;
        for (int i = 0; i < vertices[u].size; ++i)
        {
            int v = vertices[u].neighbor[i];
            vertices[v].state = 1 - vertices[v].state;
            loss_1 += 2 * vertices[v].state - 1;
        }
        for (int i = 0; i < vertices[u].size; ++i)
        {
            int v = vertices[u].neighbor[i];
            for (int j = 0; j < vertices[v].size; ++j)
            {
                int w = vertices[v].neighbor[j];
                if (w != u && !edge_matrix[u][w])
                {
                    loss_21 += 2 * (vertices[v].state ^ vertices[w].state) - 1;
                    loss_31 += (2 * (vertices[v].state ^ vertices[w].state) - 1) * (1.0 / vertices[v].size + 1.0 / vertices[w].size);
                }
            }
        }
    }

    void print()
    {
        for (vertex v : vertices)
        {
            cout << v.state << " ";
        }
        cout << endl;
    }
} g(0);

void gradient_descent()
{
    bool is_converged = false;
    while (!is_converged)
    {
        is_converged = true;
        for (int i = 0; i < len; ++i)
        {
            sum_num++;
            g.light_up(i);
            double temp_loss_num = g.loss_1;
            if (temp_loss_num < loss_num)
            {
                loss_num = temp_loss_num;
                is_converged = false;
            }
            else
            {
                g.light_up(i);
            }
        }
    }
}

void main_algorithm()
{
    g = graph(len);
    loss = 100000;
    num = 0;
    sum_num = 0;
    temp_p = p;
    double T = T_0;
    while (true)
    {
        int x = rand() % len;
        g.light_up(x);
        double temp_loss = g.loss_1 + temp_p * g.loss_31; // change this to change the loss function
        if (loss >= temp_loss)
        {
            loss = temp_loss;
            num = 0;
        }
        else if (exp((loss - temp_loss) / T) > ((double)rand() / RAND_MAX))
        {
            loss = temp_loss;
            num = 0;
        }
        else
        {
            g.light_up(x);
        }
        num++;
        sum_num++;
        T *= cooling_rate;
        if (T < T_min)
        {
            break;
        }
        if (T < 1 && T > 0.1)
        {
            temp_p *= 1 - 2e-5;
        }
        if (T < 0.1)
        {
            temp_p *= 1 - 1e-4;
        }
    }
    loss_num = g.loss_1;
    gradient_descent();
    loss = g.loss_1 + temp_p * g.loss_31; // change this to change the loss function
    // print_light();
    cout << "len: " << len << "\tp: " << p << "\tloss_num: " << loss_num << " \tloss: " << loss << " \tsum_num: " << sum_num << endl;
    sum += loss_num;
    sum2 += loss_num * loss_num;
}

void repeat_algorithm()
{
    sum = 0;
    sum2 = 0;
    is_stable = false;
    for (repeat_times = 1; repeat_times <= repeat_times_bound; ++repeat_times)
    {
        main_algorithm();
        if (repeat_times >= 10)
            if (sum2 / (repeat_times * (repeat_times - 1)) - ((sum * sum) / (repeat_times * repeat_times * (repeat_times - 1))) < (sum / repeat_times * 0.002) * (sum / repeat_times * 0.002))
            {
                is_stable = true;
                break;
            }
    }
}

void extract_last_entry(int &len_value, double &p_value)
{
    ifstream infile("light_output.txt");
    if (!infile)
    {
        cerr << "无法打开文件 light_output.txt" << endl;
        return;
    }

    string line;
    string last_line;
    while (getline(infile, line))
    {
        if (!line.empty())
        {
            last_line = line;
        }
    }

    infile.close();

    if (last_line.empty())
    {
        cerr << "文件为空或没有有效数据" << endl;
        return;
    }

    istringstream iss(last_line);
    string part;
    while (getline(iss, part, '\t'))
    {
        if (part.find("len:") != string::npos)
        {
            len_value = stod(part.substr(part.find(":") + 1));
        }
        else if (part.find("p:") != string::npos)
        {
            p_value = stod(part.substr(part.find(":") + 1));
        }
    }
}

int main()
{
    srand(time(0));

    extract_last_entry(init_len, init_p);

    ofstream outfile("light_output.txt", std::ios::app);
    if (!outfile)
    {
        cerr << "无法打开文件 light_output.txt" << endl;
        return 1;
    }

    for (len = 1000; len <= 1000; len += 10)
    {
        for (p = p_min; p <= p_max + 1e-6; p += p_step)
        {
            if (is_init)
            {
                p = init_p + p_step;
                is_init = false;
                if (p > p_max)
                    break;
            }
            repeat_algorithm();
            outfile << "len: " << len << "\tp: " << p << "\taverage loss_num: " << sum / repeat_times;
            if (is_stable)
            {
                outfile << endl;
            }
            else
            {
                outfile << "\tnot stable: " << (sqrt(sum2 / (repeat_times * (repeat_times - 1)) - ((sum * sum) / (repeat_times * repeat_times * (repeat_times - 1)))) / (sum / repeat_times)) << endl;
            }
        }
    }

    // outfile << "这是写入文件的第一行。" << endl;
    // outfile << "这是写入文件的第二行。" << endl;

    outfile.close();

    return 0;
}
