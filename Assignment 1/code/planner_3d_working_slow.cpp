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

#include <boost/functional/hash.hpp>
#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>


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

#define NUMOFDIRS 9



struct Point
{
    public:
    int x = 0;
    int y = 0;
    int t = 0;

    Point(){}

    Point(int a, int b, int c) : x(a), y(b), t(c) {}
    
    
    void set(int a, int b, int c)
    {
        this->x = a;
        this->y = b;
        this->t = c;

    }

    bool operator==(const Point& other) const
    {
        return((this->x == other.x) && (this->y == other.y) && (this->t == other.t));
    }

    bool operator!=(const Point& other)const
    {
        return((this->x != other.x) || (this->y != other.y) || (this->t != other.t));
    }

    friend std::size_t hash_value(Point const& p)
    {
        std::size_t seed = 0;
        boost::hash_combine(seed,p.x);
        boost::hash_combine(seed,p.y);
        boost::hash_combine(seed,p.t);
        return seed;
    }
};

struct Node
{
    Point coordinate;
    Point parent;
    double f = 1e9;
    double g = 1e9;
    double h = 1e9;

    Node() {}

    void set_f()
    {
        f = g + h;
    }
};

// struct PointHash {

//     size_t operator()(Point p) const {
//         std::size_t seed = 0;
//         boost::hash_combine(seed,p.x);
//         boost::hash_combine(seed,p.y);
//         boost::hash_combine(seed,p.t);
//         return seed;
// };

class Compare
{
public:
    bool operator() (Node* a, Node* b)
    {
        if (a->f != b->f) return (a->f)>(b->f);
        else return ((a->h) > (b->h));
    }
};

// Global variables declaration
std::stack<Point> global_path;

class a_star
{
    public:
    // cmp = [](Node* a, Node*b){return (a->f)>(b->f);};
    std::priority_queue<Node*,std::vector<Node*>,Compare> open;
    boost::unordered_set<Point> closed;
    int x_size, y_size;
    boost::unordered_set<std::pair<int,int>> goal_list;
    Point goalpose;
    double* target_traj;
    int target_steps;
    Point startpose;
    double* map;
    int dX[NUMOFDIRS] = {-1, -1, -1,  0, 0,  0,  1, 1, 1};
    int dY[NUMOFDIRS] = {-1,  0,  1, -1, 0,  1, -1, 0, 1};
    // double c_s[NUMOFDIRS] = {sqrt(2),1,sqrt(2),1,1,1,sqrt(2),1,sqrt(2)};
    double c_s[NUMOFDIRS] = {1,1,1,1,1,1,1,1,1};
    int collision_thresh;
    boost::unordered_map<Point,Node> my_map;

    
    a_star(int size_x, int size_y, int robotposeX, int robotposeY, double* tg_traj, int c_time, int t_steps, double* global_map, int coll_thresh) : 
    x_size(size_x), y_size(size_y), map(global_map), collision_thresh(coll_thresh), target_traj(tg_traj), target_steps(t_steps)
    {
        this->startpose = Point(robotposeX-1,robotposeY-1,0);
        goalpose = Point(target_traj[target_steps-1]-1,target_traj[2*target_steps-1]-1,target_steps);
        for (int i=0;i<t_steps;i++) //list of all goal poses
        {
            this->goal_list.insert({tg_traj[i]-1,tg_traj[i+t_steps]-1}); //cpp is 0 indexed
        }
        Node temp;
        temp.coordinate = startpose;     
        temp.g = 0;
        temp.h = euc_dist(startpose);
        temp.set_f();
        my_map[startpose] = temp;
        open.push(&my_map[startpose]);
    }

    double euc_dist(Point a)
    {
        return (double)5*sqrt((a.x - target_traj[target_steps-1])*(a.x - target_traj[target_steps-1]) + (a.y - target_traj[2*target_steps-1])*(a.y - target_traj[2*target_steps-1]) + (a.t - target_steps)*(a.t - target_steps));    }

    void compute_path()
    {
        
        while(!(open.empty()))
        {
            Node* curr = open.top();
            // mexPrintf("x = %d, y = %d, t = %d, f = %f \n",curr->coordinate.x,curr->coordinate.y,curr->coordinate.t,curr->f);
            if (curr->coordinate == goalpose)
            {
                break;
            }

            
            closed.insert(curr->coordinate);
            open.pop();
            // mexPrintf("Size of open: %d \n",open.size());
            //mexPrintf("Size of closed: %d \n",closed.size());

            for(int dir = 0; dir < NUMOFDIRS; dir++)
            {
                Point new_loc = Point(curr->coordinate.x + dX[dir],curr->coordinate.y + dY[dir],curr->coordinate.t+1);

                if (new_loc.x >= 0 && new_loc.x < x_size && new_loc.y >= 0 && new_loc.y < y_size) //Within Bounds
                {
                    if (closed.count(new_loc) == 0) //Not in closed list
                    {
                        int t_c = map[GETMAPINDEX(new_loc.x+1,new_loc.y+1,x_size,y_size)];
                        if ((t_c >= 0) && (t_c < collision_thresh))  //if free
                        {
                            Node temp;
                            temp.coordinate = new_loc;
                            temp.g = curr->g + t_c * c_s[dir]; //1 for Time
                            temp.h = euc_dist(new_loc);
                            temp.set_f();
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
            
        }
    }

    void make_path() //For now just return 1 action
    {
        std::stack<Point> path;
        Node curr = my_map[goalpose];
        path.push(curr.coordinate);
        while (curr.parent != startpose)
        {
            curr = my_map[curr.parent];
            path.push(curr.coordinate);
        }

        global_path = path;
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

    int goalposeX = (int) target_traj[curr_time];
    int goalposeY = (int) target_traj[curr_time+target_steps];
    // printf("robot: %d %d;\n", robotposeX, robotposeY);
    // printf("goal: %d %d;\n", goalposeX, goalposeY);

    a_star mystar(x_size,y_size, robotposeX, robotposeY, target_traj, curr_time, target_steps, map, collision_thresh);
    if (global_path.empty())
    {
        mystar.compute_path();
        mystar.make_path();
    }

    action_ptr[0] = global_path.top().x+1;
    action_ptr[1] = global_path.top().y+1; //matlab is 1 indexed
    global_path.pop();
    
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