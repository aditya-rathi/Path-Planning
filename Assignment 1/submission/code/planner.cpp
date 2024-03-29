/*=================================================================
 *
 * planner.c
 *
 *=================================================================*/
#include <math.h>
#include <mex.h>

#include <vector>
#include <algorithm>
#include <iterator>
#include <array>
#include <queue>
#include <stack>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <utility> 
#include <ctime>


/* Input Arguments */
#define	MAP_IN                  prhs[0]
#define	ROBOT_IN                prhs[1]
#define	TARGET_TRAJ             prhs[2]
#define	TARGET_POS              prhs[3]
#define	CURR_TIME               prhs[4]
#define	COLLISION_THRESH        prhs[5]


/* Output Arguments */
#define	ACTION_OUT              plhs[0]

//access to the map is shifted to account for 0-based indexing in the map, whereas
//1-based indexing in matlab (so, robotpose and goalpose are 1-indexed)
#define GETMAPINDEX(X, Y, XSIZE, YSIZE) ((Y-1)*XSIZE + (X-1))

#if !defined(MAX)
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define	MIN(A, B)	((A) < (B) ? (A) : (B))
#endif

#define NUMOFDIRS 8

struct Point
{
    int x = 0;
    int y = 0;

    Point(){}

    Point(int a, int b) : x(a), y(b) {}
    
    
    void set(int a, int b)
    {
        this->x = a;
        this->y = b;
    }
};

struct Node
{
    std::pair<int,int> coordinate;
    std::pair<int,int> parent;
    double f = 1e9;
    double g = 1e9;
    double h = 1e9;
    int dist = 1e9; //number of moves made

    void set_f()
    {
        f = g + h;
    }
};

struct IntPairHash {
    static_assert(sizeof(int) * 2 == sizeof(size_t));

    size_t operator()(std::pair<int,int> p) const noexcept {
        return size_t(p.first) << 32 | p.second;
    }
};

class Compare
{
public:
    bool operator() (Node* a, Node* b)
    {
        return (a->f)>(b->f);
    }
};

void copyStack(std::stack<std::pair<int,int>>& stack1, std::stack<std::pair<int,int>>& stack2) {
    stack2 = stack1;
}

//Global Vars
std::stack<std::pair<int,int>>* path;
std::stack<std::pair<int,int>>* temp_path;
int mygoal_x, mygoal_y;
double cumel_g = 0;
double best_g;

class a_star
{
    public:
    // cmp = [](Node* a, Node*b){return (a->f)>(b->f);};
    std::priority_queue<Node*,std::vector<Node*>,Compare> open;
    std::unordered_set<std::pair<int,int>,IntPairHash> closed;
    std::unordered_map<std::pair<int,int>,Node,IntPairHash> my_map;
    int x_size, y_size;
    std::pair<int,int> goalpose;
    std::pair<int,int> startpose;
    double* map;
    int dX[NUMOFDIRS] = {-1, -1, -1,  0,  0,  1, 1, 1};
    int dY[NUMOFDIRS] = {-1,  0,  1, -1,  1, -1, 0, 1};
    double c_s[NUMOFDIRS] = {1,1,1,1,1,1,1,1};
    int collision_thresh;
    int target_steps;
    int curr_time;
    bool path_success;
    clock_t start_time;
    

    
    a_star(int size_x, int size_y, int robotposeX, int robotposeY, int goalposeX, int goalposeY, double* global_map, int coll_thresh, int t_steps, int c_time, clock_t s_time) : 
    x_size(size_x), y_size(size_y), map(global_map), collision_thresh(coll_thresh), target_steps(t_steps), curr_time(c_time), start_time(s_time)
    {
        this->startpose = std::make_pair(robotposeX-1,robotposeY-1);
        this->goalpose = std::make_pair(goalposeX-1,goalposeY-1); //cpp is 0 indexed
        Node temp;
        temp.coordinate = startpose;     
        temp.g = 0;
        temp.h = euc_dist(startpose,goalpose);
        temp.set_f();
        temp.dist = 0;
        my_map[startpose] = temp;
        open.push(&my_map[startpose]);
        path_success  = false;
    }

    static bool compare_f_val(const Node* a, const Node* b)
    {
        if(a->f<b->f) return true;
        return false;
    }

    double euc_dist(std::pair<int,int> a, std::pair<int,int> b)
    {
        return (double)sqrt(((a.first-b.first)*(a.first-b.first) + (a.second-b.second)*(a.second-b.second)));
    }

    void compute_path()
    {
        
        while(!(open.empty()))
        {
            if (curr_time!=0)
            {
                if ((double(clock() - start_time)/double(CLOCKS_PER_SEC)) > 0.8) break;
            }

            Node* curr = open.top();
            if (curr->coordinate == goalpose) 
            {
                path_success = true;
                break;
            }

            closed.insert(curr->coordinate);
            open.pop();

            for(int dir = 0; dir < NUMOFDIRS; dir++) 
            {
                if(curr->dist < (target_steps - curr_time)) expand_succ(curr,dir);  
            }
        }
    }

    void expand_succ(Node* curr, int dir)
    {
        std::pair<int,int> new_loc = std::make_pair(curr->coordinate.first + dX[dir],curr->coordinate.second + dY[dir]);

        if (new_loc.first >= 0 && new_loc.first < x_size && new_loc.second >= 0 && new_loc.second < y_size) //Within Bounds
        {
            if (closed.count(new_loc) ==0) //Not in closed list
            {
                int t_c = map[GETMAPINDEX(new_loc.first+1,new_loc.second+1,x_size,y_size)];
                if ((t_c >= 0) && (t_c < collision_thresh))  //if free
                {
                    Node temp;
                    temp.coordinate = new_loc;
                    temp.g = curr->g + t_c * c_s[dir];
                    temp.h = euc_dist(new_loc,goalpose);
                    temp.set_f();
                    temp.dist = curr->dist + 1;
                    temp.parent = curr->coordinate;
                    auto it = my_map.find(new_loc);
                    if (it==my_map.end()) 
                    {
                        my_map[new_loc] = temp;
                        open.push(&my_map[new_loc]);
                        
                    }
                    else 
                    {
                        if (it->second.f > temp.f)
                        {
                            it->second.f = temp.f; 
                            it->second.g = temp.g;
                            it->second.parent = temp.parent;
                        }
                    }
                
                }
                else //If not free add to closed list
                {
                    closed.insert(new_loc);
                }
            }
        }
    }

    void make_path() //For now just return 1 action
    {
        temp_path = new std::stack<std::pair<int,int>>;
        Node curr = my_map[goalpose];
        best_g = curr.g + cumel_g + (target_steps - curr_time - curr.dist);
        while (curr.parent != startpose)
        {
            temp_path->push(curr.coordinate);
            curr = my_map[curr.parent];
        }
        temp_path->push(curr.coordinate);
        path = temp_path;
    }

    
};

static void planner(
        double*	map,
        int collision_thresh,
        int x_size,
        int y_size,
        int robotposeX,
        int robotposeY,
        int target_steps,
        double* target_traj,
        int targetposeX,
        int targetposeY,
        int curr_time,
        double* action_ptr
        )
{
    clock_t start_time = clock();
    // 8-connected grid
    int dX[NUMOFDIRS] = {-1, -1, -1,  0,  0,  1, 1, 1};
    int dY[NUMOFDIRS] = {-1,  0,  1, -1,  1, -1, 0, 1};
    
    int goalposeX = (int) target_traj[target_steps-1];
    int goalposeY = (int) target_traj[target_steps-1+target_steps];

    // mexPrintf("goal x: %d, y = %d \n", mygoal_x,mygoal_y);
    // mexPrintf("robot x: %d, y = %d \n",robotposeX,robotposeY);
    
    // Get current target position
    if (curr_time==0)
    {
        mygoal_x = goalposeX;
        mygoal_y = goalposeY;
        if(robotposeX != goalposeX || robotposeY != goalposeY)
        {
            a_star mystar(x_size, y_size, robotposeX, robotposeY, goalposeX, goalposeY, map, collision_thresh,target_steps,curr_time,start_time);
            mystar.compute_path();
            mystar.make_path();
        }
        else
        {
            action_ptr[0] = goalposeX;
            action_ptr[1] = goalposeY;
            return;
        }
    }
    else
    {
        //Backtrack the target path
        if (target_steps - curr_time > curr_time)
        {
            int goalposeX = (int) target_traj[target_steps - 1 - curr_time];
            int goalposeY = (int) target_traj[target_steps - 1 + target_steps - curr_time];
            
            
            a_star mystar(x_size,y_size, robotposeX, robotposeY, goalposeX, goalposeY, map, collision_thresh,target_steps,curr_time,start_time);
            mystar.compute_path();
            if ((mystar.my_map[{goalposeX,goalposeY}].g + cumel_g + (target_steps - curr_time - mystar.my_map[{goalposeX,goalposeY}].dist)) < best_g) 
            {
                // mexPrintf("%f \n",mystar.my_map[{goalposeX,goalposeY}].g + cumel_g);
                if (mystar.path_success)
                {
                mystar.make_path();
                mygoal_x = goalposeX;
                mygoal_y = goalposeY;
                }
            }
        }    
    }
    if(robotposeX != mygoal_x || robotposeY != mygoal_y)
    {
        std::pair<int,int> act = path->top();
        cumel_g+=map[GETMAPINDEX(act.first+1,act.second+1,x_size,y_size)];
        action_ptr[0] = act.first+1;
        action_ptr[1] = act.second+1; //matlab is 1 indexed

        path->pop();
    }
    else
    {
        // delete(path);
        action_ptr[0] = mygoal_x;
        action_ptr[1] = mygoal_y;
    }
    // mexPrintf("%f \n",double(clock() - start_time)/double(CLOCKS_PER_SEC));
    return;
}

// prhs contains input parameters (4):
// 1st is matrix with all the obstacles
// 2nd is a row vector <x,y> for the robot position
// 3rd is a matrix with the target trajectory
// 4th is an integer C, the collision threshold for the map
// plhs should contain output parameters (1):
// 1st is a row vector <dx,dy> which corresponds to the action that the robot should make
void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[] )
        
{
    
    /* Check for proper number of arguments */
    if (nrhs != 6) {
        mexErrMsgIdAndTxt( "MATLAB:planner:invalidNumInputs",
                "Six input arguments required.");
    } else if (nlhs != 1) {
        mexErrMsgIdAndTxt( "MATLAB:planner:maxlhs",
                "One output argument required.");
    }
    
    /* get the dimensions of the map and the map matrix itself*/
    int x_size = mxGetM(MAP_IN);
    int y_size = mxGetN(MAP_IN);
    double* map = mxGetPr(MAP_IN);
    
    /* get the dimensions of the robotpose and the robotpose itself*/
    int robotpose_M = mxGetM(ROBOT_IN);
    int robotpose_N = mxGetN(ROBOT_IN);
    if(robotpose_M != 1 || robotpose_N != 2){
        mexErrMsgIdAndTxt( "MATLAB:planner:invalidrobotpose",
                "robotpose vector should be 1 by 2.");
    }
    double* robotposeV = mxGetPr(ROBOT_IN);
    int robotposeX = (int)robotposeV[0];
    int robotposeY = (int)robotposeV[1];
    
    /* get the dimensions of the goalpose and the goalpose itself*/
    int targettraj_M = mxGetM(TARGET_TRAJ);
    int targettraj_N = mxGetN(TARGET_TRAJ);
    
    if(targettraj_M < 1 || targettraj_N != 2)
    {
        mexErrMsgIdAndTxt( "MATLAB:planner:invalidtargettraj",
                "targettraj vector should be M by 2.");
    }
    double* targettrajV = mxGetPr(TARGET_TRAJ);
    int target_steps = targettraj_M;
    
    /* get the current position of the target*/
    int targetpose_M = mxGetM(TARGET_POS);
    int targetpose_N = mxGetN(TARGET_POS);
    if(targetpose_M != 1 || targetpose_N != 2){
        mexErrMsgIdAndTxt( "MATLAB:planner:invalidtargetpose",
                "targetpose vector should be 1 by 2.");
    }
    double* targetposeV = mxGetPr(TARGET_POS);
    int targetposeX = (int)targetposeV[0];
    int targetposeY = (int)targetposeV[1];
    
    /* get the current timestep the target is at*/
    int curr_time = mxGetScalar(CURR_TIME);
    
    /* Create a matrix for the return action */ 
    ACTION_OUT = mxCreateNumericMatrix( (mwSize)1, (mwSize)2, mxDOUBLE_CLASS, mxREAL); 
    double* action_ptr = (double*) mxGetData(ACTION_OUT);
    
    /* Get collision threshold for problem */
    int collision_thresh = (int) mxGetScalar(COLLISION_THRESH);
    
    /* Do the actual planning in a subroutine */
    planner(map, collision_thresh, x_size, y_size, robotposeX, robotposeY, target_steps, targettrajV, targetposeX, targetposeY, curr_time, &action_ptr[0]);
    // printf("DONE PLANNING!\n");
    return;   
}