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
    Point coordinate;
    Node* parent;
    double f = 1e9;
    double g = 1e9;
    double h = 1e9;

    Node(){}

    Node(int a, int b) : coordinate(a,b) {}

    void set_f()
    {
        f = g + h;
    }
};

// bool operator==(const Node* lhs, const Node* rhs)
// {
//     if (lhs == rhs) return true;
//     return false;
// }

class a_star
{
    public:
    std::vector<Node*> open;
    std::vector<Node*> closed;
    int x_size, y_size;
    Node goalpose;
    Node startpose;
    double* g_vals;
    double* map;
    int dX[NUMOFDIRS] = {-1, -1, -1,  0,  0,  1, 1, 1};
    int dY[NUMOFDIRS] = {-1,  0,  1, -1,  1, -1, 0, 1};
    int collision_thresh;
    
    a_star(int size_x, int size_y, int robotposeX, int robotposeY, int goalposeX, int goalposeY, double* global_map, int coll_thresh) : 
    x_size(size_x), y_size(size_y), map(global_map), collision_thresh(coll_thresh)
    {
        this->goalpose.coordinate.set(goalposeX, goalposeY);
        this->startpose.coordinate.set(robotposeX,robotposeY);
        this->startpose.g = 0;
        this->startpose.h = euc_dist(startpose,goalpose);
        this->startpose.set_f();
        open.push_back(&startpose);
    }

    static bool compare_f_val(const Node* a, const Node* b)
    {
        if(a->f<b->f) return true;
        return false;
    }

    double euc_dist(Node a, Node b)
    {
        return (double)sqrt(((a.coordinate.x-b.coordinate.x)*(a.coordinate.x-b.coordinate.x) + (a.coordinate.y-b.coordinate.y)*(a.coordinate.y-b.coordinate.y)));
    }

    void compute_path()
    {
        std::cout<<((std::find(open.begin(),open.end(),&(this->goalpose))==open.end()) && !(open.empty()));
        while((std::find(open.begin(),open.end(),&(this->goalpose))==open.end()) && !(open.empty()))
        {
            std::sort(open.begin(), open.end(),[](Node* a, Node* b){return (a->f<b->f);}); //compare_f_val);
            Node* curr = open[0];
            mexPrintf("x: %d, y= %d \n",(*curr).coordinate.x,(*curr).coordinate.y);
            closed.push_back(curr);
            std::vector<Node*> succ;
            for(int dir = 0; dir < NUMOFDIRS; dir++)
            {
                Node* new_loc = new Node(curr->coordinate.x + dX[dir],curr->coordinate.y + dY[dir]);

                if (new_loc->coordinate.x >= 1 && new_loc->coordinate.x <= x_size && new_loc->coordinate.y >= 1 && new_loc->coordinate.y <= y_size) //Within Bounds
                {
                    if (find_if(closed.begin(),closed.end(),[new_loc](Node* a){return ((new_loc->coordinate.x == a->coordinate.x)&&(new_loc->coordinate.y == a->coordinate.y));})==closed.end())
                    {
                        if (((int)map[GETMAPINDEX(new_loc->coordinate.x,new_loc->coordinate.y,x_size,y_size)] >= 0) && ((int)map[GETMAPINDEX(new_loc->coordinate.x,new_loc->coordinate.y,x_size,y_size)] < collision_thresh))  //if free
                        {
                            if (new_loc->g > curr->g + euc_dist(*curr, *new_loc))
                            {
                                new_loc->g = curr->g + euc_dist(*curr, *new_loc);
                                new_loc->h = euc_dist(*new_loc,goalpose);
                                new_loc->set_f();
                                new_loc->parent = curr;
                                std::vector<Node*>::iterator it = std::find(open.begin(),open.end(),new_loc);
                                if (it==open.end()) open.push_back(new_loc);
                                else {(**it).f = new_loc->f; (**it).parent = new_loc->parent;}
                            }
                        }
                        else //If not free add to closed list
                        {
                            closed.push_back(new_loc);
                        }
                    }
                }
            }
            open.erase(open.begin());
        }
    }

    Node make_path() //For now just return 1 action
    {
        std::stack<Node*> path;
        Node* curr = &goalpose;
        while (curr->parent != &startpose)
        {
            path.push(curr);
            curr = curr->parent;
        }

        return *(path.top());
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

    // 8-connected grid
    int dX[NUMOFDIRS] = {-1, -1, -1,  0,  0,  1, 1, 1};
    int dY[NUMOFDIRS] = {-1,  0,  1, -1,  1, -1, 0, 1};
    
    // Get current target position

    int goalposeX = (int) target_traj[curr_time-1];
    int goalposeY = (int) target_traj[curr_time-1+target_steps];
    // printf("robot: %d %d;\n", robotposeX, robotposeY);
    // printf("goal: %d %d;\n", goalposeX, goalposeY);

    a_star mystar(x_size,y_size, robotposeX, robotposeY, goalposeX, goalposeY, map, collision_thresh);
    mystar.compute_path();
    Node act = mystar.make_path();

    action_ptr[0] = act.coordinate.x;
    action_ptr[1] = act.coordinate.y;
    
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