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
#include <limits>

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
            // cout << endl;
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
	Node* parent = NULL;
};

bool myequalDoubleArrays(double* v1, double *v2, int size) {
    for (int i = 0; i < size; ++i) {
        if (abs(v1[i]-v2[i]) > 1e-1) {
            // cout << endl;
            return false;
        }
    }
    return true;
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

int nearest_nhb(const std::vector<double*> node_list, double* curr, int numofDOFs)
{
	double best_dist = euc_dist(node_list[0],curr,numofDOFs);
	int nhb_ind = 0;

	for (int i=1;i<node_list.size();i++)
	{
		double temp_dist = euc_dist(node_list[i],curr,numofDOFs);
		if (temp_dist<best_dist)
		{
			best_dist = temp_dist;
			nhb_ind = i;
		}
	}
	return nhb_ind;
}

std::pair<bool,double*> add_elem(const std::vector<double*> node_list, int nearest_nhb_ind, double* curr, int numofDOFs, double eps_start, double* map, int x_size, int y_size)
{
	double eps = eps_start;
	double* nhb_coord = node_list[nearest_nhb_ind];
	double dist = euc_dist(nhb_coord,curr,numofDOFs);
	double unit_vec[5];
	double* new_node = (double*) malloc(numofDOFs*sizeof(double));
	
	if(dist<eps) eps = dist;

	for(int i=0;i<numofDOFs;i++) unit_vec[i] = (curr[i] - nhb_coord[i])/dist;

	for(int i=0;i<numofDOFs;i++) new_node[i] = nhb_coord[i] + unit_vec[i]*eps;

	while(!IsValidArmConfiguration(new_node,numofDOFs,map,x_size,y_size) && eps>0)
	{
		eps-=1e-2;
		for(int i=0;i<numofDOFs;i++) new_node[i] = nhb_coord[i] + unit_vec[i]*eps;
	}

	if(eps<0) return std::make_pair(false,new_node);
	return std::make_pair(true,new_node);
}

int uf_find(std::vector<int>& uf_array, int x)
{
	if (uf_array[x] != x) uf_array[x] = uf_find(uf_array,uf_array[x]);
	return uf_array[x];
}

void uf_union(std::vector<int>& uf_array, std::vector<int>& uf_rank, int src, int dst)
{
	int x = uf_find(uf_array,src);
	int y = uf_find(uf_array,dst);
	if (uf_rank[x] < uf_rank[y])
	{
		uf_array[x] = y;
	}
	else if (uf_rank[x] > uf_rank[y])
	{
		uf_array[y] = x;
	}
	else
	{
		uf_array[y] = x;
		uf_rank[x]++;
	}
}

std::vector<int> k_nearest_nhb(const std::vector<double*>& samples, double* curr, int k, int numofDOFs)
{
	std::priority_queue<std::pair<double,int>,std::vector<std::pair<double,int>>,std::greater<std::pair<double,int>>> pq;
	for (int i=0; i<samples.size(); i++)
	{
		double dist = euc_dist(samples[i],curr,numofDOFs);
		if (dist>1e-3)	pq.push({dist,i});
	}
	std::vector<int> result;
	for (int i= 0;i<k;i++)
	{
		result.push_back(pq.top().second);
		pq.pop();
	}
	pq = std::priority_queue<std::pair<double,int>,std::vector<std::pair<double,int>>,std::greater<std::pair<double,int>>>();
	return result;
}

void rrt(double* map, int x_size, int y_size, double* armstart_anglesV_rad, double* armgoal_anglesV_rad, int numofDOFs, double*** plan, int* planlength)
{
	int nodes_added = 0;
	double eps = 0.2;
	//no plan by default
	*plan = NULL;
	*planlength = 0;

	//Initialize
	std::vector<double*> node_list;
	std::vector<int> parent_list;
	node_list.push_back(armstart_anglesV_rad);
	parent_list.push_back(-1);


	//Generate random sample
	double* nhb_coord = (double*) malloc(numofDOFs*sizeof(double));

	while(!myequalDoubleArrays(nhb_coord,armgoal_anglesV_rad,numofDOFs) && nodes_added<15000)
	{
		double* curr = random_config(map, x_size, y_size, numofDOFs);
		int nearest_nhb_ind = nearest_nhb(node_list,curr,numofDOFs);
		
		std::pair<bool, double*> data;
		data = add_elem(node_list, nearest_nhb_ind, curr, numofDOFs, eps, map, x_size, y_size);
		
		if(data.first)
		{
			node_list.push_back(data.second);
			parent_list.push_back(nearest_nhb_ind);
			nodes_added++;
			nhb_coord = data.second;
		}	
	}
	
	// Reconstruct path
	int parent_ind = nearest_nhb(node_list,armgoal_anglesV_rad,numofDOFs);
	std::stack<double*> path;
	path.push(armgoal_anglesV_rad);
	
	while(parent_ind!=-1)
	{
		path.push(node_list[parent_ind]);
		parent_ind = parent_list[parent_ind];
	}

	*planlength = path.size();
	*plan = (double**) malloc((path.size())*sizeof(double*));

    for (int i = 0; i < *planlength; i++)
	{
        (*plan)[i] = (double*) malloc(numofDOFs*sizeof(double)); 
        for(int j = 0; j < numofDOFs; j++)
		{
            (*plan)[i][j] = path.top()[j];
        }
		path.pop();
	}
	
}

void rrt_connect(double* map, int x_size, int y_size, double* armstart_anglesV_rad, double* armgoal_anglesV_rad, int numofDOFs, double*** plan, int* planlength)
{
	int nodes_added = 0;
	double eps = 0.2;
	//no plan by default
	*plan = NULL;
	*planlength = 0;

	//Initialize
	std::vector<double*> fwd_node_list;
	std::vector<int> fwd_parent_list;
	std::vector<double*> back_node_list;
	std::vector<int> back_parent_list;
	bool which_tree = true; //true = fwd_tree || false = back_tree
	bool trees_connected = false;
	
	fwd_node_list.push_back(armstart_anglesV_rad);
	fwd_parent_list.push_back(-1);
	back_node_list.push_back(armgoal_anglesV_rad);
	back_parent_list.push_back(-1);


	//Generate random sample
	double* last_node_added;
	std::pair<bool, double*> data; 

	while(!trees_connected && nodes_added<15000)
	{
		// std::cout<<nodes_added<<std::endl;
		double* sample = random_config(map, x_size, y_size, numofDOFs);
		int nearest_nhb_ind;
		if(which_tree)	nearest_nhb_ind = nearest_nhb(fwd_node_list,sample,numofDOFs);
		else nearest_nhb_ind = nearest_nhb(back_node_list,sample,numofDOFs);
		
		
		if (which_tree)	data = add_elem(fwd_node_list, nearest_nhb_ind, sample, numofDOFs, eps, map, x_size, y_size);
		else data = add_elem(back_node_list, nearest_nhb_ind, sample, numofDOFs, eps, map, x_size, y_size);
		//expand current main tree
		if(data.first)
		{
			last_node_added = data.second;
			nodes_added++;

			if (which_tree)
			{	
				fwd_node_list.push_back(data.second);
				fwd_parent_list.push_back(nearest_nhb_ind);
			}
			else
			{
				back_node_list.push_back(data.second);
				back_parent_list.push_back(nearest_nhb_ind);
			}

			// Try connecting other tree
			data = std::make_pair(true,static_cast<double*>(NULL));
			while(data.first)
			{
				if (which_tree)
				{
					sample = fwd_node_list[fwd_node_list.size()-1];
					nearest_nhb_ind = nearest_nhb(back_node_list,sample,numofDOFs);
					data = add_elem(back_node_list, nearest_nhb_ind, sample, numofDOFs, eps, map, x_size, y_size);
					back_node_list.push_back(data.second);
					back_parent_list.push_back(nearest_nhb_ind);
				}
				else
				{
					sample = back_node_list[back_node_list.size()-1];
					nearest_nhb_ind = nearest_nhb(fwd_node_list,sample,numofDOFs);
					data = add_elem(fwd_node_list, nearest_nhb_ind, sample, numofDOFs, eps, map, x_size, y_size);
					fwd_node_list.push_back(data.second);
					fwd_parent_list.push_back(nearest_nhb_ind);
				}
				if(equalDoubleArrays(data.second,last_node_added,numofDOFs))
				{
					trees_connected = true;
					break;	
				} 
			}
		}	
		which_tree = !which_tree; //flip trees
	}
	
	// Reconstruct path
	std::deque<double*> path;
	int fwd_ind, back_ind;
	if(!trees_connected)
	{
		fwd_ind = nearest_nhb(fwd_node_list,armgoal_anglesV_rad,numofDOFs);
		back_ind = nearest_nhb(back_node_list,fwd_node_list[fwd_ind],numofDOFs);
	}
	else
	{
		fwd_ind = fwd_node_list.size()-1;
		back_ind = back_node_list.size()-1;
	}
	
	
	
	// Make fwd path first
	while(fwd_ind!=-1)
	{
		path.push_back(fwd_node_list[fwd_ind]);
		fwd_ind = fwd_parent_list[fwd_ind];
	}
	// Next make back path
	while(back_ind!=-1)
	{
		path.push_front(back_node_list[back_ind]);
		back_ind = back_parent_list[back_ind];
	}

	*planlength = path.size();
	*plan = (double**) malloc((*planlength)*sizeof(double*));

    for (int i = 0; i < *planlength; i++)
	{
        (*plan)[i] = (double*) malloc(numofDOFs*sizeof(double)); 
        for(int j = 0; j < numofDOFs; j++)
		{
            (*plan)[i][j] = path.back()[j];
        }
		path.pop_back();
	}
	
}

void prm(double* map, int x_size, int y_size, double* armstart_anglesV_rad, double* armgoal_anglesV_rad, int numofDOFs, double*** plan, int* planlength)
{
	int nodes_added = 0;
	double eps = 0.3;

	//no plan by default
	*plan = NULL;
	*planlength = 0;

	//Preprocessing Phase
	const int numsamples = 10000;
	std::vector<double*> samples;
	std::vector<std::vector<int>> edges(numsamples+2,std::vector<int>(numsamples+2,0));
	
	// int edges[numsamples+2][numsamples+2]; //Adjacency matrix of edges. Add 2 for start and goal
	std::vector<int> uf_array; //to store components
	std::vector<int> uf_rank(numsamples,0);

	//Generate the random samples
	for(int i=0;i<numsamples;i++) 
	{
		samples.push_back(random_config(map,x_size,y_size,numofDOFs));
		uf_array.push_back(i);
	}

	// Assuming k=3 for now
	int k = 3;

	for(int i=0;i<numsamples;i++)
	{
		double* curr = samples[i];
		if (i==4451)
		{
			std::cout<<"hello"<<std::endl;
		}

		// Create edge map		
		std::vector<int> temp = k_nearest_nhb(samples,curr,k,numofDOFs);		
		for (int j = 0;j<temp.size();j++)
		{
			
			int nhb_ind = temp[j];
			if (edges[i][nhb_ind] == 0) 
			{
				edges[i][nhb_ind] = 1;
				// Since our graph is undirected
				edges[nhb_ind][i] = 1;
			}

			// Create connected components
			uf_union(uf_array,uf_rank,i,nhb_ind); //set parent of nhb_ind to i			
		}
	}
		
	//Planning phase
	//Augment graph with start and goal
	samples.push_back(armstart_anglesV_rad); //numsamples+1
	samples.push_back(armgoal_anglesV_rad); //numsamples+2
	
	k = 5; //Check 5 closest neighbors when connecting start and goal
	std::vector<int> start_nhb = k_nearest_nhb(samples, armstart_anglesV_rad, k, numofDOFs);
	std::vector<int> goal_nhb = k_nearest_nhb(samples, armgoal_anglesV_rad, k, numofDOFs);
	std::vector<int> comps;
	for (int i=0; i<k; i++)
	{
		int start_comp = uf_find(uf_array,start_nhb[i]);
		comps.push_back(start_comp);
	}
	for (int i=0; i<k; i++)
	{
		bool connected = false;
		int goal_comp = uf_find(uf_array,goal_nhb[i]);
		for (int j = 0;j<k;j++)
		{
			if(goal_comp == comps[j])
			{
				edges[numsamples][start_nhb[j]] = 1;
				edges[start_nhb[j]][numsamples] = 1;
				edges[numsamples+1][goal_nhb[i]] = 1;
				edges[goal_nhb[i]][numsamples+1] = 1;
				connected = true;
				break;
			}
		}
		if(connected) break;
	}
	

	//Djikstra's algo for path finding
	std::vector<int> distances(numsamples+2,std::numeric_limits<int>::max());
	std::vector<bool> visited(numsamples+2,false);
	std::vector<int> parents;
	parents.reserve(numsamples+2);

	//Initialize
	distances[numsamples] = 0;
	parents[numsamples] = -1;

	//Find path
	for(int i=0;i<numsamples+2;i++)
	{
		int curr = -1;
		int curr_dist = std::numeric_limits<int>::max();
		for (int j = 0; j<numsamples+2;j++)
		{
			if(!visited[j] && distances[j]<curr_dist)
			{
				curr = j;
				curr_dist = distances[j];
			}
		}
		visited[curr] = true;
		for (int j = 0; j<numsamples+2;j++)
		{
			int edgecost = edges[curr][j];
			if (edgecost>0 && (curr_dist+edgecost < distances[j]))
			{
				distances[j] = curr_dist + edgecost;
				parents[j] = curr;
			}
		}
		if (distances[numsamples+1]<std::numeric_limits<int>::max()) break; //Goal has been reached
	}

	//Make path
	std::stack<double*> path;
	double* curr = samples[numsamples+1];
	int parent = parents[numsamples+1];
	while (parent!=-1)
	{
		path.push(curr);
		curr = samples[parent];
		parent = parents[parent];
	}
	path.push(curr);
	*planlength = path.size();
	*plan = (double**) malloc((path.size())*sizeof(double*));

    for (int i = 0; i < *planlength; i++)
	{
        (*plan)[i] = (double*) malloc(numofDOFs*sizeof(double)); 
        for(int j = 0; j < numofDOFs; j++)
		{
            (*plan)[i][j] = path.top()[j];
        }
		path.pop();
	}
}



static void planner(
			double* map,
			int x_size,
			int y_size,
			double* armstart_anglesV_rad,
			double* armgoal_anglesV_rad,
            int numofDOFs,
            double*** plan,
            int* planlength,
			int whichPlanner)
{
	//seed the engine
	eng.seed(16782);

	switch (whichPlanner)
	{
	case 0:
		std::cout<<"Using rrt"<<std::endl;
		rrt(map,x_size,y_size,armstart_anglesV_rad,armgoal_anglesV_rad,numofDOFs,plan,planlength);
		break;
	
	case 1:
		std::cout<<"Using rrt-connect"<<std::endl;
		rrt_connect(map,x_size,y_size,armstart_anglesV_rad,armgoal_anglesV_rad,numofDOFs,plan,planlength);
		break;

	case 3:
		std::cout<<"Using prm"<<std::endl;
		prm(map,x_size,y_size,armstart_anglesV_rad,armgoal_anglesV_rad,numofDOFs,plan,planlength);
		break;

	default:
		throw runtime_error("Invalid planner selected!\n");
		break;
	}
	std::cout<<"Path found successfully"<<std::endl;
	//    
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
	planner(map, x_size, y_size, startPos, goalPos, numOfDOFs, &plan, &planlength, whichPlanner);

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
