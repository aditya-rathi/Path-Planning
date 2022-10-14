/*=================================================================
 *
 * planner.c
 *
 *=================================================================*/
#include <math.h>
#include <random>
#include <vector>
#include <array>
#include <algorithm>

#include <tuple>
#include <string>
#include <stdexcept>
#include <regex> // For regex and split logic
#include <iostream> // cout, endl
#include <fstream> // For reading/writing files
#include <assert.h> 

#include <random>
#include <queue>

/* Input Arguments */
#define	MAP_IN      prhs[0]
#define	ARMSTART_IN	prhs[1]
#define	ARMGOAL_IN     prhs[2]
#define	PLANNER_ID_IN     prhs[3]

/* Planner Ids */
#define RRT         0
#define RRTCONNECT  1
#define RRTSTAR     2
#define PRM         3

/* Output Arguments */
#define	PLAN_OUT	plhs[0]
#define	PLANLENGTH_OUT	plhs[1]

#define GETMAPINDEX(X, Y, XSIZE, YSIZE) (Y*XSIZE + X)

#if !defined(MAX)
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define	MIN(A, B)	((A) < (B) ? (A) : (B))
#endif

#define PI 3.141592654

//the length of each link in the arm
#define LINKLENGTH_CELLS 10

// Some potentially helpful imports
using std::vector;
using std::array;
using std::string;
using std::runtime_error;
using std::tuple;
using std::make_tuple;
using std::tie;
using std::cout;
using std::endl;

/// @brief 
/// @param filepath 
/// @return map, x_size, y_size
tuple<double*, int, int> loadMap(string filepath) {
	std::FILE *f = fopen(filepath.c_str(), "r");
	if (f) {
	}
	else {
		printf("Opening file failed! \n");
		throw runtime_error("Opening map file failed!");
	}
	int height, width;
	if (fscanf(f, "height %d\nwidth %d\n", &height, &width) != 2) {
		throw runtime_error("Invalid loadMap parsing map metadata");
	}
	
	////// Go through file and add to m_occupancy
	double* map = new double[height*width];

	double cx, cy, cz;
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			char c;
			do {
				if (fscanf(f, "%c", &c) != 1) {
					throw runtime_error("Invalid parsing individual map data");
				}
			} while (isspace(c));
			if (!(c == '0')) { 
				map[y+x*width] = 1; // Note transposed from visual
			} else {
				map[y+x*width] = 0;
			}
		}
	}
	fclose(f);
	return make_tuple(map, width, height);
}

// Splits string based on deliminator
vector<string> split(const string& str, const string& delim) {   
		// https://stackoverflow.com/questions/14265581/parse-split-a-string-in-c-using-string-delimiter-standard-c/64886763#64886763
		const std::regex ws_re(delim);
		return { std::sregex_token_iterator(str.begin(), str.end(), ws_re, -1), std::sregex_token_iterator() };
}


double* doubleArrayFromString(string str) {
	vector<string> vals = split(str, ",");
	double* ans = new double[vals.size()];
	for (int i = 0; i < vals.size(); ++i) {
		ans[i] = std::stod(vals[i]);
	}
	return ans;
}

bool equalDoubleArrays(double* v1, double *v2, int size) {
    for (int i = 0; i < size; ++i) {
        if (abs(v1[i]-v2[i]) > 1e-3) {
            cout << endl;
            return false;
        }
    }
    return true;
}

typedef struct {
	int X1, Y1;
	int X2, Y2;
	int Increment;
	int UsingYIndex;
	int DeltaX, DeltaY;
	int DTerm;
	int IncrE, IncrNE;
	int XIndex, YIndex;
	int Flipped;
} bresenham_param_t;


void ContXY2Cell(double x, double y, short unsigned int* pX, short unsigned int *pY, int x_size, int y_size) {
	double cellsize = 1.0;
	//take the nearest cell
	*pX = (int)(x/(double)(cellsize));
	if( x < 0) *pX = 0;
	if( *pX >= x_size) *pX = x_size-1;

	*pY = (int)(y/(double)(cellsize));
	if( y < 0) *pY = 0;
	if( *pY >= y_size) *pY = y_size-1;
}


void get_bresenham_parameters(int p1x, int p1y, int p2x, int p2y, bresenham_param_t *params) {
	params->UsingYIndex = 0;

	if (fabs((double)(p2y-p1y)/(double)(p2x-p1x)) > 1)
		(params->UsingYIndex)++;

	if (params->UsingYIndex)
		{
			params->Y1=p1x;
			params->X1=p1y;
			params->Y2=p2x;
			params->X2=p2y;
		}
	else
		{
			params->X1=p1x;
			params->Y1=p1y;
			params->X2=p2x;
			params->Y2=p2y;
		}

	 if ((p2x - p1x) * (p2y - p1y) < 0)
		{
			params->Flipped = 1;
			params->Y1 = -params->Y1;
			params->Y2 = -params->Y2;
		}
	else
		params->Flipped = 0;

	if (params->X2 > params->X1)
		params->Increment = 1;
	else
		params->Increment = -1;

	params->DeltaX=params->X2-params->X1;
	params->DeltaY=params->Y2-params->Y1;

	params->IncrE=2*params->DeltaY*params->Increment;
	params->IncrNE=2*(params->DeltaY-params->DeltaX)*params->Increment;
	params->DTerm=(2*params->DeltaY-params->DeltaX)*params->Increment;

	params->XIndex = params->X1;
	params->YIndex = params->Y1;
}

void get_current_point(bresenham_param_t *params, int *x, int *y) {
	if (params->UsingYIndex) {
        *y = params->XIndex;
        *x = params->YIndex;
        if (params->Flipped)
            *x = -*x;
    }
	else {
        *x = params->XIndex;
        *y = params->YIndex;
        if (params->Flipped)
            *y = -*y;
    }
}

int get_next_point(bresenham_param_t *params) {
	if (params->XIndex == params->X2) {
        return 0;
    }
	params->XIndex += params->Increment;
	if (params->DTerm < 0 || (params->Increment < 0 && params->DTerm <= 0))
		params->DTerm += params->IncrE;
	else {
        params->DTerm += params->IncrNE;
        params->YIndex += params->Increment;
	}
	return 1;
}



int IsValidLineSegment(double x0, double y0, double x1, double y1, double*	map,
			 int x_size, int y_size) {
	bresenham_param_t params;
	int nX, nY; 
	short unsigned int nX0, nY0, nX1, nY1;

	//printf("checking link <%f %f> to <%f %f>\n", x0,y0,x1,y1);
		
	//make sure the line segment is inside the environment
	if(x0 < 0 || x0 >= x_size ||
		x1 < 0 || x1 >= x_size ||
		y0 < 0 || y0 >= y_size ||
		y1 < 0 || y1 >= y_size)
		return 0;

	ContXY2Cell(x0, y0, &nX0, &nY0, x_size, y_size);
	ContXY2Cell(x1, y1, &nX1, &nY1, x_size, y_size);

	//printf("checking link <%d %d> to <%d %d>\n", nX0,nY0,nX1,nY1);

	//iterate through the points on the segment
	get_bresenham_parameters(nX0, nY0, nX1, nY1, &params);
	do {
		get_current_point(&params, &nX, &nY);
		if(map[GETMAPINDEX(nX,nY,x_size,y_size)] == 1)
			return 0;
	} while (get_next_point(&params));

	return 1;
}

int IsValidArmConfiguration(double* angles, int numofDOFs, double*	map,
			 int x_size, int y_size) {
    double x0,y0,x1,y1;
    int i;
		
	 //iterate through all the links starting with the base
	x1 = ((double)x_size)/2.0;
	y1 = 0;
	for(i = 0; i < numofDOFs; i++){
		//compute the corresponding line segment
		x0 = x1;
		y0 = y1;
		x1 = x0 + LINKLENGTH_CELLS*cos(2*PI-angles[i]);
		y1 = y0 - LINKLENGTH_CELLS*sin(2*PI-angles[i]);

		//check the validity of the corresponding line segment
		if(!IsValidLineSegment(x0,y0,x1,y1,map,x_size,y_size))
			return 0;
	}    
	return 1;
}

//Global random num engine
std::mt19937 eng;

struct Node
{
	vector<double> coord;
	Node* left = NULL;
	Node* right = NULL;
	Node* parent = NULL;
};

double euc_dist(vector<double> a, vector<double> b)
{
	if(a.size()!=b.size()) return -1;
	double result = 0;
	for(int i=0;i<a.size();i++) result+=((a[i]-b[i])*(a[i]-b[i]));
	return sqrt(result);
}

double euc_dist(double* a, double* b, int size)
{
	double result=0;
	for(int i=0;i<size;i++) result+=((a[i]-b[i])*(a[i]-b[i]));
	return sqrt(result);	
}

double* random_config(double* map, int x_size, int y_size, int numofDOFs)
{
	std::uniform_real_distribution<double> num_gen(0,2*PI);
	double* config = (double*) malloc(numofDOFs*sizeof(double));
	for (int i = 0; i<numofDOFs; i++)
	{
		config[i] = num_gen(eng); 
	}
	while(!IsValidArmConfiguration(config,numofDOFs,map,x_size,y_size))
	{
		for (int i = 0; i<numofDOFs; i++)
		{
			config[i] = num_gen(eng); 
		}
	}
	return config;
}

Node* nearest_nhb(Node* head, double* curr_arr, int numofDOFs, double* result)
{
	vector<double> curr = vector<double>(curr_arr,curr_arr+numofDOFs);
	int axis = 0;
	int depth = 0;
	vector<double> nearest_nhb = head->coord;
	double best_dist = euc_dist(head->coord,curr);
	Node* parent = head;

	while(head)
	{
		double temp_dist = euc_dist(head->coord,curr);
		if (temp_dist<best_dist)
		{
			best_dist = temp_dist;
			nearest_nhb = head->coord;
			parent = head;
		}

		if (curr[axis]<head->coord[axis]) head = head->left;
		else if (curr[axis]>head->coord[axis]) head = head->right;
		else break;

		depth++;
		axis = depth%numofDOFs;
	}
	std::copy(nearest_nhb.begin(),nearest_nhb.end(),result);
	return parent; 
}

void add_elem(Node* head, Node* temp, double* nhb_coord, double* curr, int numofDOFs, double eps_start, double* map, int x_size, int y_size)
{
	double eps = eps_start;
	double dist = euc_dist(nhb_coord,curr,numofDOFs);
	double unit_vec[5];
	for(int i=0;i<numofDOFs;i++)
	{
		unit_vec[i] = (curr[i] - nhb_coord[i])/dist;
	}

	for(int i=0;i<numofDOFs;i++)
	{
		temp->coord[i] = nhb_coord[i] + unit_vec[i]*eps;
	}

	while(!IsValidArmConfiguration(&temp->coord[0],numofDOFs,map,x_size,y_size) && eps>0)
	{
		eps-=1e-2;
		for(int i=0;i<numofDOFs;i++)
		{
			temp->coord[i] = nhb_coord[i] + unit_vec[i]*eps;
		}
	}

	if(eps<0) return;
	
	int axis = 0;
	int depth = 0;
	Node* parent = head;
	while(head)
	{
		parent = head;
		if(temp->coord[axis]<head->coord[axis]) head = head->left;
		else head = head->right;
		depth++;
		axis = depth%numofDOFs;
	}

	axis = (depth-1)%numofDOFs;
	if(temp->coord[axis] < parent->coord[axis]) parent->left = temp;
	else parent->right = temp;
	return;
}

void rrt(double* map, int x_size, int y_size, double* armstart_anglesV_rad, double* armgoal_anglesV_rad, int numofDOFs, double*** plan, int* planlength)
{
	int nodes_added = 0;
	double eps = 0.5;
	//no plan by default
	*plan = NULL;
	*planlength = 0;

	//Initialize K-d tree
	Node* head = new Node();
	head->coord = std::vector<double>(armstart_anglesV_rad,armstart_anglesV_rad+numofDOFs);


	//Generate random sample
	double* nhb_coord = (double*) malloc(numofDOFs*sizeof(double));
	for(int i=0;i<numofDOFs; nhb_coord[i] = head->coord[i], i++);

	Node* last_added;

	while(!equalDoubleArrays(nhb_coord,armgoal_anglesV_rad,numofDOFs))
	{
		double* curr = random_config(map, x_size, y_size, numofDOFs);
		Node* parent = nearest_nhb(head,curr,numofDOFs,nhb_coord);
		last_added = new Node();
		last_added->coord.assign(5,0);
		last_added->parent = parent;

		add_elem(head, last_added, nhb_coord,curr,numofDOFs,eps, map, x_size, y_size);
		nodes_added++;
	}
	
	// Reconstruct path
	std::queue<double*> path;
	path.push(armgoal_anglesV_rad);
	while(last_added!=head)
	{
		last_added = last_added->parent;
		path.push(&last_added->coord[0]);
	}
	*plan = (double**) malloc((path.size()+1)*sizeof(double*));
	for(int j=0; j<numofDOFs; j++) (*plan)[0][j] = armstart_anglesV_rad[j];
    for (int i = 0; i < (path.size()+1); i++){
        (*plan)[i] = (double*) malloc(numofDOFs*sizeof(double)); 
        for(int j = 0; j < numofDOFs; j++){
            (*plan)[i][j] = *(path.front()+j);
        }
		path.pop();
	}
	*planlength = path.size()+1;
}



static void planner(
			double* map,
			int x_size,
			int y_size,
			double* armstart_anglesV_rad,
			double* armgoal_anglesV_rad,
            int numofDOFs,
            double*** plan,
            int* planlength)
{
	//seed the engine
	eng.seed(16782);

	rrt(map,x_size,y_size,armstart_anglesV_rad,armgoal_anglesV_rad,numofDOFs,plan,planlength);

	//no plan by default
	// *plan = NULL;
	// *planlength = 0;
		
    // //for now just do straight interpolation between start and goal checking for the validity of samples

    // double distance = 0;
    // int i,j;
    // for (j = 0; j < numofDOFs; j++){
    //     if(distance < fabs(armstart_anglesV_rad[j] - armgoal_anglesV_rad[j]))
    //         distance = fabs(armstart_anglesV_rad[j] - armgoal_anglesV_rad[j]);
    // }
    // int numofsamples = (int)(distance/(PI/20));
    // if(numofsamples < 2){
    //     printf("The arm is already at the goal\n");
    //     return;
    // }
	// int countNumInvalid = 0;
    // *plan = (double**) malloc(numofsamples*sizeof(double*));
    // for (i = 0; i < numofsamples; i++){
    //     (*plan)[i] = (double*) malloc(numofDOFs*sizeof(double)); 
    //     for(j = 0; j < numofDOFs; j++){
    //         (*plan)[i][j] = armstart_anglesV_rad[j] + ((double)(i)/(numofsamples-1))*(armgoal_anglesV_rad[j] - armstart_anglesV_rad[j]);
    //     }
    //     if(!IsValidArmConfiguration((*plan)[i], numofDOFs, map, x_size, y_size)) {
	// 		++countNumInvalid;
    //     }
    // }
	// printf("Linear interpolation collided at %d instances across the path\n", countNumInvalid);
    // *planlength = numofsamples;
    
    return;
}


/** Your final solution will be graded by an grading script which will
 * send the default 6 arguments:
 *    map, numOfDOFs, commaSeparatedStartPos, commaSeparatedGoalPos, 
 *    whichPlanner, outputFilePath
 * An example run after compiling and getting the planner.out executable
 * >> ./planner.out map1.txt 5 1.57,0.78,1.57,0.78,1.57 0.392,2.35,3.14,2.82,4.71 0 output.txt
 * See the hw handout for full information.
 * If you modify this for testing (e.g. to try out different hyper-parameters),
 * make sure it can run with the original 6 commands.
 * Programs that do not will automatically get a 0.
 * */
int main(int argc, char** argv) {
	double* map;
	int x_size, y_size;

	tie(map, x_size, y_size) = loadMap(argv[1]);
	const int numOfDOFs = std::stoi(argv[2]);
	double* startPos = doubleArrayFromString(argv[3]);
	double* goalPos = doubleArrayFromString(argv[4]);
	int whichPlanner = std::stoi(argv[5]);
	string outputFile = argv[6];

	if(!IsValidArmConfiguration(startPos, numOfDOFs, map, x_size, y_size)||
			!IsValidArmConfiguration(goalPos, numOfDOFs, map, x_size, y_size)) {
		throw runtime_error("Invalid start or goal configuration!\n");
	}

	///////////////////////////////////////
	//// Feel free to modify anything below. Be careful modifying anything above.

	double** plan = NULL;
	int planlength = 0;
	planner(map, x_size, y_size, startPos, goalPos, numOfDOFs, &plan, &planlength);

	//// Feel free to modify anything above.
	//// If you modify something below, please change it back afterwards as my 
	//// grading script will not work and you will recieve a 0.
	///////////////////////////////////////

    // Your solution's path should start with startPos and end with goalPos
    if (!equalDoubleArrays(plan[0], startPos, numOfDOFs) || 
    	!equalDoubleArrays(plan[planlength-1], goalPos, numOfDOFs)) {
		throw std::runtime_error("Start or goal position not matching");
	}

	/** Saves the solution to output file
	 * Do not modify the output log file output format as it is required for visualization
	 * and for grading.
	 */
	std::ofstream m_log_fstream;
	m_log_fstream.open(outputFile, std::ios::trunc); // Creates new or replaces existing file
	if (!m_log_fstream.is_open()) {
		throw std::runtime_error("Cannot open file");
	}
	m_log_fstream << argv[1] << endl; // Write out map name first
	/// Then write out all the joint angles in the plan sequentially
	for (int i = 0; i < planlength; ++i) {
		for (int k = 0; k < numOfDOFs; ++k) {
			m_log_fstream << plan[i][k] << ",";
		}
		m_log_fstream << endl;
	}
}
