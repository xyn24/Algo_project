#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <string>
#include <regex>

using namespace std;

const double p_max = 3;
int len = 50;
const int repeat_times_bound = 1000;
int repeat_times = 0;
double sum = 0;
double sum2 = 0;
int len_x = len;
int len_y = len;
const double T_0 = 100;
const double cooling_rate = 1 - 1e-4;
const double T_min = 0.01;
const int d_xy[5][2] = {{0, 0}, {0, 1}, {0, -1}, {1, 0}, {-1, 0}};
int light_stage[200][200];
double loss = 100000;
double loss_num = 100000;
int num = 0;
int sum_num = 0;
double p = 0, temp_p = 0;
int init_len = 20;
double init_p = 0;
bool is_init = false;  //change this to initialize the p
bool is_stable = false;

void light_up_kernel(int light_stage[200][200], int x_, int y_)
{
    for (int i = 0; i < 5; ++i)
    {
        int new_x = x_ + d_xy[i][0];
        int new_y = y_ + d_xy[i][1];
        if (new_x >= 0 && new_x < len_x && new_y >= 0 && new_y < len_y)
        {
            light_stage[new_x][new_y] = 1 - light_stage[new_x][new_y];
        }
    }
}

void light_up(int x_, int y_)
{
    light_up_kernel(light_stage, x_, y_);
}

void generate_random_light()
{
    for (int i = 0; i < len_x; ++i)
    {
        for (int j = 0; j < len_y; ++j)
        {
            light_stage[i][j] = rand() % 2;
        }
    }
}

double find_new_loss_1()
{
    double temp_loss = 0;
    for (int i = 0; i < len_x; ++i)
    {
        for (int j = 0; j < len_y; ++j)
        {
            temp_loss += light_stage[i][j];
        }
    }
    return temp_loss;
}

double find_new_loss_1(int new_p)
{
    double temp_loss = 0;
    for (int i = 0; i < len_x; ++i)
    {
        for (int j = 0; j < len_y; ++j)
        {
            temp_loss += light_stage[i][j];
        }
    }
    return temp_loss;
}

double find_new_loss_2(int new_p)
{
    double temp_loss = 0;
    for (int i = 0; i < len_x - 1; ++i)
    {
        for (int j = 0; j < len_y; ++j)
        {
            temp_loss += light_stage[i][j] * light_stage[i + 1][j];
        }
    }
    for (int i = 0; i < len_x; ++i)
    {
        for (int j = 0; j < len_y - 1; ++j)
        {
            temp_loss += light_stage[i][j] * light_stage[i][j + 1];
        }
    }
    return -new_p * temp_loss + find_new_loss_1();
}

double find_new_loss_3(int new_p)
{
    double temp_loss = 0;
    for (int i = 0; i < len_x - 1; ++i)
    {
        for (int j = 0; j < len_y; ++j)
        {
            temp_loss += light_stage[i][j]^light_stage[i + 1][j];
        }
    }
    for (int i = 0; i < len_x; ++i)
    {
        for (int j = 0; j < len_y - 1; ++j)
        {
            temp_loss += light_stage[i][j]^light_stage[i][j + 1];
        }
    }
    return new_p * temp_loss + find_new_loss_1();
}

void print_light()
{
    for (int i = 0; i < len_y; ++i)
    {
        for (int j = 0; j < len_x; ++j)
        {
            cout << light_stage[i][j] << " ";
        }
        cout << endl;
    }
}

void get_light_input()
{
    for (int i = 0; i < len_y; ++i)
    {
        for (int j = 0; j < len_x; ++j)
        {
            cin >> light_stage[i][j];
        }
    }
}

void gradient_descent()
{
    bool is_converged = false;
    while (!is_converged)
    {
        is_converged = true;
        for (int i = 0; i < len_x; ++i) 
            for (int j = 0; j < len_y; ++j)
            {
                sum_num++;
                light_up(i, j);
                double temp_loss_num = find_new_loss_1();
                light_up(i, j);
                if (temp_loss_num < loss_num)
                {
                    light_up(i, j);
                    loss_num = temp_loss_num;
                    is_converged = false;
                }
            }
    }
}

void main_algorithm()
{
    generate_random_light();
    loss = 100000;
    num = 0;
    sum_num = 0;
    temp_p = p;
    double T = T_0;
    while (true)
    {
        int x = rand() % len_x;
        int y = rand() % len_y;
        light_up(x, y);
        double temp_loss = find_new_loss_3(temp_p); // change this to change the loss function
        light_up(x, y);
        if (loss >= temp_loss)
        {
            light_up(x, y);
            loss = temp_loss;
            num = 0;
        }
        else if (exp((loss - temp_loss) / T) > ((double)rand() / RAND_MAX))
        {
            light_up(x, y);
            loss = temp_loss;
            num = 0;
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
            temp_p *= 1-2e-5;
        }
        if (T < 0.1)
        {
            temp_p *= 1-1e-4;
        }
    }
    loss_num = find_new_loss_1();
    gradient_descent();
    loss = find_new_loss_3(temp_p); // change this to change the loss function
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
            if (sum2 / (repeat_times*(repeat_times-1)) - ((sum*sum) / (repeat_times*repeat_times*(repeat_times-1))) < (sum / repeat_times * 0.005) * (sum / repeat_times * 0.005))
            {
                is_stable = true;
                break;
            }
    }
}

int main()
{
    srand(time(0));

    ifstream infile("light_output.txt");
    if (!infile) {
        cerr << "无法打开文件 light_output.txt" << endl;
        return 1;
    }

    string line;
    regex pattern(R"(len: (\d+)\s+p: (\d+\.\d+))");
    smatch matches;

    while (getline(infile, line)) {
        if (regex_search(line, matches, pattern)) {
            init_len = stoi(matches[1].str());
            init_p = double(stod(matches[2].str()))+0.05;
        }
    }

    infile.close();

    ofstream outfile("light_output.txt", std::ios::app);
    if (!outfile) {
        cerr << "无法打开文件 light_output.txt" << endl;
        return 1;
    }
    
    for (len = 50; len<=50; len+=10)
    {
        len_x = len;
        len_y = len;
        for (p = p_max; p <= p_max+1e-6; p += 0.05)
        {
            if (is_init)
            {
                p = init_p;
                is_init = false;
                if (p>p_max)
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
                outfile << "\tnot stable: " << (sqrt(sum2 / (repeat_times*(repeat_times-1)) - ((sum*sum) / (repeat_times*repeat_times*(repeat_times-1))))/(sum / repeat_times)) << endl;
            }
        }
    }

    // outfile << "这是写入文件的第一行。" << endl;
    // outfile << "这是写入文件的第二行。" << endl;

    outfile.close();

    return 0;
}
