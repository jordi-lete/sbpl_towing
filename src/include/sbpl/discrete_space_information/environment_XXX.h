/*
 * Copyright (c) 2008, Maxim Likhachev
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the Carnegie Mellon University nor the names of its
 *       contributors may be used to endorse or promote products derived from
 *       this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef __ENVIRONMENT_XXXLAT_H_
#define __ENVIRONMENT_XXXLAT_H_

#include <cstdio>
#include <vector>
#include <sstream>
#include <unordered_map>

#include <sbpl/discrete_space_information/environment.h>
#include <sbpl/utils/utils.h>

// Define to test against in client code. Signals that Set2DBlockSize and
// Set2DBucketSize are available in EnvironmentXXXLATTICE
#define SBPL_CUSTOM_2D_OPTIONS 1

//eight-connected grid
#define XXXLAT_DXYWIDTH 8
#define EnvXXXLAT_DEFAULTOBSTHRESH 254	//see explanation of the value below
//maximum number of states for storing them into lookup (as opposed to hash)
#define SBPL_XXXLAT_MAXSTATESFORLOOKUP 100000000
//definition of theta orientations
//0 - is aligned with X-axis in the positive direction (1,0 in polar coordinates)
//theta increases as we go counterclockwise
//number of theta values - should be power of 2
#define XXXLAT_THETADIRS 16
//number of actions per x,y,theta state
//decrease, increase, same angle while moving plus decrease, increase angle while standing.
#define XXXLAT_DEFAULT_ACTIONWIDTH 5
#define XXXLAT_COSTMULT_MTOMM 1000

class CMDPSTATE;
class MDPConfig;
class SBPL2DGridSearch;

struct TrailerState {
    double x;
    double y;
    double theta1;
    double theta2; // Used only for four wheel trailers
    TrailerState() : x(0), y(0), theta1(0), theta2(0) {}
    TrailerState(double x_, double y_, double theta1_, double theta2_ = 0.0) : x(x_), y(y_), theta1(theta1_), theta2(theta2_) {}
};

struct PathState {
    double x;
    double y;
    double theta;
    double time;
};

struct PathHistory {
    std::vector<PathState> positions;
    TrailerState currentTrailer;
};

struct EnvXXXLATAction_t
{
    unsigned char aind; //index of the action (unique for given starttheta)
    char starttheta;
    char dX;
    char dY;
    char endtheta;
    double time;
    unsigned int cost;
    std::vector<sbpl_2Dcell_t> intersectingcellsV;
    //start at 0,0,starttheta and end at endcell in continuous domain with half-bin less to account for 0,0 start
    std::vector<sbpl_xy_theta_pt_t> intermptV;
    //start at 0,0,starttheta and end at endcell in discrete domain
    std::vector<sbpl_xy_theta_cell_t> interm3DcellsV;

 int motprimID;
 double turning_radius;

};

struct EnvXXXLATHashEntry_t
{
    int stateID;
    int X;
    int Y;
    char Theta;
    int iteration;
};

struct SBPL_xxx_mprimitive
{
    int motprimID;
    unsigned char starttheta_c;
    int additionalactioncostmult;
    sbpl_xy_theta_cell_t endcell;
    double turning_radius;
    //intermptV start at 0,0,starttheta and end at endcell in continuous
    //domain with half-bin less to account for 0,0 start
    std::vector<sbpl_xy_theta_pt_t> intermptV;
};

//variables that dynamically change (e.g., array of states, ...)
struct EnvironmentXXXLAT_t
{
    int startstateid;
    int goalstateid;

    bool bInitialized;

    //any additional variables
};

//configuration parameters
struct EnvXXXLATConfig_t
{
    int EnvWidth_c;
    int EnvHeight_c;
    int NumThetaDirs;
    int StartX_c;
    int StartY_c;
    int StartTheta;
    double StartTheta1;
    double StartTheta2;
    int EndX_c;
    int EndY_c;
    int EndTheta;
    double goaltol_x;
    double goaltol_y;
    double goaltol_theta;
    unsigned char** Grid2D;

    double R0;
    double F1;
    double F2;
    int num_pivots;

    std::vector<double> ThetaDirs;
    double StartTheta_rad;
    double EndTheta_rad;
    double min_turning_radius_m;

    // the value at which and above which cells are obstacles in the maps sent from outside
    // the default is defined above
    unsigned char obsthresh;

    // the value at which and above which until obsthresh (not including it)
    // cells have the nearest obstacle at distance smaller than or equal to
    // the inner circle of the robot. In other words, the robot is definitely
    // colliding with the obstacle, independently of its orientation
    // if no such cost is known, then it should be set to obsthresh (if center
    // of the robot collides with obstacle, then the whole robot collides with
    // it independently of its rotation)
    unsigned char cost_inscribed_thresh;

    // the value at which and above which until cost_inscribed_thresh (not including it) cells
    // **may** have a nearest osbtacle within the distance that is in between
    // the robot inner circle and the robot outer circle
    // any cost below this value means that the robot will NOT collide with any
    // obstacle, independently of its orientation
    // if no such cost is known, then it should be set to 0 or -1 (then no cell
    // cost will be lower than it, and therefore the robot's footprint will
    // always be checked)
    int cost_possibly_circumscribed_thresh; // it has to be integer, because -1 means that it is not provided.

    double nominalvel_mpersecs;

    //double nominalangvel_radpersecs;

    double timetoturn45degsinplace_secs;

    double cellsize_m;

    int dXY[XXXLAT_DXYWIDTH][2];

    //array of actions, ActionsV[i][j] - jth action for sourcetheta = i
    EnvXXXLATAction_t** ActionsV;
    //PredActionsV[i] - vector of pointers to the actions that result in a state with theta = i
    std::vector<EnvXXXLATAction_t*>* PredActionsV;

    int actionwidth; //number of motion primitives
    std::vector<SBPL_xxx_mprimitive> mprimV;

    std::vector<sbpl_2Dpt_t> FootprintPolygon;
    std::vector<sbpl_2Dpt_t> TrailerPolygon;
};

class EnvXXXLAT_InitParms
{
public:
    unsigned int numThetas;
    const unsigned char* mapdata;
    double startx;
    double starty;
    double starttheta;
    double starttheta1;
    double starttheta2;
    double goalx;
    double goaly;
    double goaltheta;
    double goaltol_x;
    double goaltol_y;
    double goaltol_theta;
    double R0;
    double F1;
    double F2;
    int num_pivots;
};

/** \brief 3D (x,y,theta) planning using lattice-based graph problem. For
 *         general structure see comments on parent class DiscreteSpaceInformation
 *         For info on lattice-based planning used here, you can check out the paper:
 *         Maxim Likhachev and Dave Ferguson, " Planning Long Dynamically-Feasible
 *         Maneuvers for Autonomous Vehicles", IJRR'09
 */
class EnvironmentXXXLATTICE : public DiscreteSpaceInformation
{
public:
    EnvironmentXXXLATTICE();

    virtual double calculateContAngleDiff(double startangle, double endangle) = 0;

    /**
     * \brief initialization of environment from file. See .cfg files for
     *        examples it also takes the perimeter of the robot with respect to some
     *        reference point centered at x=0,y=0 and orientation = 0 (along x axis).
     *        The perimeter is defined in meters as a sequence of vertices of a
     *        polygon defining the perimeter. If vector is of zero size, then robot
     *        is assumed to be point robot (you may want to inflate all obstacles by
     *        its actual radius) Motion primitives file defines the motion primitives
     *        available to the robot
     */
    virtual bool InitializeEnv(const char* sEnvFile, const std::vector<sbpl_2Dpt_t>& perimeterptsV,
                               const std::vector<sbpl_2Dpt_t>& trailer_perimeterptsV, const char* sMotPrimFile);

    //  ??  virtual bool InitializeEnv(const ENVNAVXYTHETAVELOLAT_InitParms & params);

    /**
     * \brief see comments on the same function in the parent class
     */
    virtual bool InitializeEnv(const char* sEnvFile);

    /**
     * \brief way to set up various parameters. For a list of parameters, see
     *        the body of the function - it is pretty straightforward
     */
    bool SetEnvParameter(const char* parameter, int value);

    /**
     * \brief returns the value of specific parameter - see function body for the list of parameters
     */
    int GetEnvParameter(const char* parameter);

    /**
     * \brief see comments on the same function in the parent class
     */
    bool InitializeMDPCfg(MDPConfig *MDPCfg);

    /**
     * \brief see comments on the same function in the parent class
     */
    virtual int GetFromToHeuristic(int FromStateID, int ToStateID) = 0;

    /**
     * \brief see comments on the same function in the parent class
     */
    virtual int GetGoalHeuristic(int stateID) = 0;

    /**
     * \brief see comments on the same function in the parent class
     */
    virtual int GetStartHeuristic(int stateID) = 0;

    /**
     * \brief see comments on the same function in the parent class
     */
    virtual void SetAllActionsandAllOutcomes(CMDPSTATE* state) = 0;

    /**
     * \brief see comments on the same function in the parent class
     */
    void SetAllPreds(CMDPSTATE* state);

    /**
     * \brief see comments on the same function in the parent class
     */
    virtual void GetSuccs(int SourceStateID, std::vector<int>* SuccIDV, std::vector<int>* CostV);

    /**
     * \brief see comments on the same function in the parent class
     */
    virtual void GetPreds(int TargetStateID, std::vector<int>* PredIDV, std::vector<int>* CostV) = 0;

    /**
     * \brief see comments on the same function in the parent class
     */
    void EnsureHeuristicsUpdated(bool bGoalHeuristics);

    /**
     * \brief see comments on the same function in the parent class
     */
    void PrintEnv_Config(FILE* fOut);

    /**
     * \brief set the block size for the 2D heuristic. A block size of 1 is the default and will result in
     *        a single cell in the 2D heuristic search corresponding to 1 cell from the source map.
     *        A block size of 2 will result in a single cell in the 2D heuristic search corresponding to
     *        a 2x2 set of blocks from the source map with a cost of the max of the 2x2 source cells.
     */
    void Set2DBlockSize(int BlockSize);

    /**
     * @brief Set2DBucketSize Set the initial size of the CSlidingBuckets used for the fringe priority list
     * @param BucketSize
     */
    void Set2DBucketSize(int BucketSize);

    double DiscTheta2ContNew(int theta) const;

    int ContTheta2DiscNew(double theta) const;

    double DiscTheta2ContFromSet(int theta) const;

    int ContTheta2DiscFromSet(double theta) const;

    int normalizeDiscAngle(int theta) const;

    /**
     * \brief initialize environment. Gridworld is defined as matrix A of size width by height.
     *        So, internally, it is accessed as A[x][y] with x ranging from 0 to width-1 and and y from 0 to height-1
     *        Each element in A[x][y] is unsigned char. A[x][y] = 0 corresponds to
     *        fully traversable and cost is just Euclidean distance
     *        The cost of transition between two neighboring cells is
     *        EuclideanDistance*(max(A[sourcex][sourcey],A[targetx][targety])+1)
     *        f A[x][y] >= obsthresh, then in the above equation it is assumed to be infinite.
     *        The cost also incorporates the length of a motion primitive and its cost_multiplier (see getcost function)
     *        mapdata is a pointer to the values of A. If it is null, then A is
     *        initialized to all zeros. Mapping is: A[x][y] = mapdata[x+y*width]
     *        start/goal are given by startx, starty, starttheta, goalx,goaly, goaltheta in meters/radians.
     *        If they are not known yet, just set them to 0. Later setgoal/setstart can be executed
     *        finally obsthresh defined obstacle threshold, as mentioned above
     *        goaltolerances are currently ignored
     *        for explanation of perimeter, see comments for InitializeEnv function that reads all from file
     *        cellsize is discretization in meters
     *        nominalvel_mpersecs is assumed velocity of vehicle while moving forward in m/sec
     *        timetoturn45degsinplace_secs is rotational velocity in secs/45 degrees turn
     */
    virtual bool InitializeEnv(int width, int height,
                               /** if mapdata is NULL the grid is initialized to all freespace */
                               const unsigned char* mapdata,
                               double startx, double starty, double starttheta, double starttheta1, double starttheta2,
                               double goalx, double goaly, double goaltheta,
                               double goaltol_x, double goaltol_y, double goaltol_theta,
                               const std::vector<sbpl_2Dpt_t>& perimeterptsV, const std::vector<sbpl_2Dpt_t>& trailer_perimeterptsV,
                               double cellsize_m, double nominalvel_mpersecs, double timetoturn45degsinplace_secs,
                               unsigned char obsthresh, const char* sMotPrimFile,
                               double R0, double F1, double F2, int num_pivots);

    /**
     * \brief Same as the above InitializeEnv except that only the parameters
     *        that area really needed are required.  The additional (optional)
     *        parameters may be given in the params object (including the ability to
     *        specify the number of thetas)
     */
    virtual bool InitializeEnv(int width, int height, const std::vector<sbpl_2Dpt_t> & perimeterptsV, const std::vector<sbpl_2Dpt_t> & trailer_perimeterptsV,
                               double cellsize_m, double nominalvel_mpersecs, double timetoturn45degsinplace_secs,
                               unsigned char obsthresh, const char* sMotPrimFile, EnvXXXLAT_InitParms params);

    /**
     * \brief update the traversability of a cell<x,y>
     */
    bool UpdateCost(int x, int y, unsigned char newcost);

    /**
     * \brief re-setting the whole 2D map
     *        transform from linear array mapdata to the 2D matrix used internally: Grid2D[x][y] = mapdata[x+y*width]
     */
    bool SetMap(const unsigned char* mapdata);

    /**
     * \brief this function fill in Predecessor/Successor states of edges whose costs changed
     *        It takes in an array of cells whose traversability changed, and
     *        returns (in vector preds_of_changededgesIDV) the IDs of all
     *        states that have outgoing edges that go through the changed
     *        cells
     */
    virtual void GetPredsofChangedEdges(std::vector<nav2dcell_t> const * changedcellsV,
                                        std::vector<int> *preds_of_changededgesIDV) = 0;
    /**
     * \brief same as GetPredsofChangedEdges, but returns successor states.
     *        Both functions need to be present for incremental search
     */
    virtual void GetSuccsofChangedEdges(std::vector<nav2dcell_t> const * changedcellsV,
                                        std::vector<int> *succs_of_changededgesIDV) = 0;

    /**
     * returns true if cell is untraversable
     */
    bool IsObstacle(int x, int y);

    /**
     * \brief returns false if robot intersects obstacles or lies outside of
     *        the map. Note this is pretty expensive operation since it computes the
     *        footprint of the robot based on its x,y,theta
     */
    bool IsValidConfiguration(int X, int Y, int Theta);

    /**
     * \brief returns environment parameters. Useful for creating a copy environment
     */
    void GetEnvParms(int *size_x, int *size_y, double* startx, double* starty, double* starttheta,
                             double* goalx, double* goaly, double* goaltheta, double* cellsize_m,
                             double* nominalvel_mpersecs, double* timetoturn45degsinplace_secs,
                             unsigned char* obsthresh, std::vector<SBPL_xxx_mprimitive>* motionprimitiveV);

    /**
     * \brief returns environment parameters. Useful for creating a copy environment
     */
    void GetEnvParms(int *size_x, int *size_y, int* num_thetas, double* startx, double* starty,
                             double* starttheta, double* goalx, double* goaly, double* goaltheta, double* cellsize_m,
                             double* nominalvel_mpersecs, double* timetoturn45degsinplace_secs,
                             unsigned char* obsthresh, std::vector<SBPL_xxx_mprimitive>* motionprimitiveV);

    /**
     * \brief get internal configuration data structure
     */
    const EnvXXXLATConfig_t* GetEnvNavConfig();

    ~EnvironmentXXXLATTICE();

    /**
     * \brief prints time statistics
     */
    void PrintTimeStat(FILE* fOut);

    /**
     * \brief returns the cost corresponding to the cell <x,y>
     */
    unsigned char GetMapCost(int x, int y);

    /**
     * \brief returns true if cell is within map
     */
    bool IsWithinMapCell(int X, int Y);

    /**
     * \brief Transform a pose into discretized form. The angle 'pth' is
     *        considered to be valid if it lies between -2pi and 2pi (some
     *        people will prefer 0<=pth<2pi, others -pi<pth<=pi, so this
     *        compromise should suit everyone).
     *
     * \note Even if this method returns false, you can still use the
     *       computed indices, for example to figure out how big your map
     *       should have been.
     *
     * \return true if the resulting indices lie within the grid bounds
     *         and the angle was valid.
     */
    bool PoseContToDisc(double px, double py, double pth, int &ix, int &iy, int &ith) const;

    /** \brief Transform grid indices into a continuous pose. The computed
     *         angle lies within 0<=pth<2pi.
     *
     * \note Even if this method returns false, you can still use the
     *      computed indices, for example to figure out poses that lie
     *      outside of your current map.
     *
     * \return true if all the indices are within grid bounds.
     */
    bool PoseDiscToCont(int ix, int iy, int ith, double &px, double &py, double &pth) const;

    /**
     * \brief prints environment variables for debugging
     */
    virtual void PrintVars() { }

protected:
    int GetActionCost(int SourceX, int SourceY, int SourceTheta, EnvXXXLATAction_t* action);

    //member data
    EnvXXXLATConfig_t EnvXXXCfg;
    EnvironmentXXXLAT_t EnvXXXLAT;
    std::vector<sbpl_xy_theta_cell_t> affectedsuccstatesV; //arrays of states whose outgoing actions cross cell 0,0
    std::vector<sbpl_xy_theta_cell_t> affectedpredstatesV; //arrays of states whose incoming actions cross cell 0,0
    int iteration;
    int blocksize; // 2D block size
    int bucketsize; // 2D bucket size
    bool bUseNonUniformAngles;

    // const double R0;
    // const double F1;
    // const double F2;

    std::unordered_map<int, PathHistory> stateToPathHistory;

    //2D search for heuristic computations
    bool bNeedtoRecomputeStartHeuristics; //set whenever grid2Dsearchfromstart needs to be re-executed
    bool bNeedtoRecomputeGoalHeuristics; //set whenever grid2Dsearchfromgoal needs to be re-executed
    SBPL2DGridSearch* grid2Dsearchfromstart; //computes h-values that estimate distances from start x,y to all cells
    SBPL2DGridSearch* grid2Dsearchfromgoal; //computes h-values that estimate distances to goal x,y from all cells

    void ReadConfiguration(FILE* fCfg);

    void InitializeEnvConfig(std::vector<SBPL_xxx_mprimitive>* motionprimitiveV);

    bool CheckQuant(FILE* fOut);

    void SetConfiguration(int width, int height,
                                  /** if mapdata is NULL the grid is initialized to all freespace */
                                  const unsigned char* mapdata,
                                  int startx, int starty, int starttheta, double starttheta1, double starttheta2,
                                  int goalx, int goaly, int goaltheta,
                                  double cellsize_m, double nominalvel_mpersecs, double timetoturn45degsinplace_secs,
                                  const std::vector<sbpl_2Dpt_t> & robot_perimeterV,
                                  const std::vector<sbpl_2Dpt_t> & trailer_perimeterV);

    bool InitGeneral(std::vector<SBPL_xxx_mprimitive>* motionprimitiveV);
    void PrecomputeActionswithBaseMotionPrimitive(std::vector<SBPL_xxx_mprimitive>* motionprimitiveV);
    void PrecomputeActionswithCompleteMotionPrimitive(std::vector<SBPL_xxx_mprimitive>* motionprimitiveV);

    virtual void InitializeEnvironment() = 0;

    void ComputeHeuristicValues();

    bool IsValidCell(int X, int Y);

    void CalculateFootprintForPose(sbpl_xy_theta_pt_t pose, std::vector<sbpl_2Dcell_t>* footprint);
    void CalculateFootprintForPose(sbpl_xy_theta_pt_t pose, std::vector<sbpl_2Dcell_t>* footprint,
                                           const std::vector<sbpl_2Dpt_t>& FootprintPolygon);
    void RemoveSourceFootprint(sbpl_xy_theta_pt_t sourcepose, std::vector<sbpl_2Dcell_t>* footprint);
    void RemoveSourceFootprint(sbpl_xy_theta_pt_t sourcepose, std::vector<sbpl_2Dcell_t>* footprint,
                                       const std::vector<sbpl_2Dpt_t>& FootprintPolygon);

    virtual void GetSuccs(int SourceStateID, std::vector<int>* SuccIDV, std::vector<int>* CostV,
                          std::vector<EnvXXXLATAction_t*>* actionindV = NULL) = 0;
    virtual int GetTrueCost(int parentID, int childID) = 0;
    virtual bool isGoal(int id) = 0;

    double EuclideanDistance_m(int X1, int Y1, int X2, int Y2);

    void ComputeReplanningData();
    void ComputeReplanningDataforAction(EnvXXXLATAction_t* action);

    bool ReadMotionPrimitives(FILE* fMotPrims);
    bool ReadinMotionPrimitive(SBPL_xxx_mprimitive* pMotPrim, FILE* fIn);
    bool ReadinCell(sbpl_xy_theta_cell_t* cell, FILE* fIn);
    bool ReadinPose(sbpl_xy_theta_pt_t* pose, FILE* fIn);

    void PrintHeuristicValues();

    bool calculateTrailerFromPath(const std::vector<PathState>& path, TrailerState& finalTrailer);
    bool calculateTrailerTransition(double startX, double startY, double startTheta, double endX, double endY, double endTheta,
                                            double time, const TrailerState& startTrailer, TrailerState& endTrailer);
    bool IsValidTrailerConfiguration(double x, double y, double theta, int& cost);
};

class EnvironmentXXXLAT : public EnvironmentXXXLATTICE
{
public:
    EnvironmentXXXLAT()
    {
        HashTableSize = 0;
        Coord2StateIDHashTable = NULL;
        Coord2StateIDHashTable_lookup = NULL;
    }

    ~EnvironmentXXXLAT();

    double calculateContAngleDiff(double startangle, double endangle);

    /**
     * \brief sets start in meters/radians
     */
    int SetStart(double x, double y, double theta);
    void SetTrailerStart(double theta1_rad, double theta2_rad);

    /**
     * \brief sets goal in meters/radians
     */
    int SetGoal(double x, double y, double theta);

    /**
     * \brief sets goal tolerance. (Note goal tolerance is ignored currently)
     */
    void SetGoalTolerance(double tol_x, double tol_y, double tol_theta) { /**< not used yet */ }

    /**
     * \brief returns state coordinates of state with ID=stateID
     */
    void GetCoordFromState(int stateID, int& x, int& y, int& theta) const;

    /**
     * \brief returns stateID for a state with coords x,y,theta
     */
    int GetStateFromCoord(int x, int y, int theta);

    /**
     * \brief returns the actions / motion primitives of the passed path.
     */
    void GetActionsFromStateIDPath(std::vector<int>* stateIDPath,
                                           std::vector<EnvXXXLATAction_t>* action_list);

    /** \brief converts a path given by stateIDs into a sequence of
     *         coordinates. Note that since motion primitives are short actions
     *         represented as a sequence of points,
     *         the path returned by this function contains much more points than the
     *         number of points in the input path. The returned coordinates are in
     *         meters,meters,radians
     */
    void ConvertStateIDPathintoXYThetaPath(std::vector<int>* stateIDPath,
                                                   std::vector<sbpl_xy_theta_pt_t>* xythetaPath, std::vector<sbpl_xy_theta_pt_t>* trailerPath);

    /**
     * \brief prints state info (coordinates) into file
     */
    void PrintState(int stateID, bool bVerbose, FILE* fOut = NULL);

    /**
     * \brief returns all predecessors states and corresponding costs of actions
     */
    void GetPreds(int TargetStateID, std::vector<int>* PredIDV, std::vector<int>* CostV);

    /**
     * \brief returns all successors states, costs of corresponding actions
     *        and pointers to corresponding actions, each of which is a motion
     *        primitive
     *        if actionindV is NULL, then pointers to actions are not returned
     */
    void GetSuccs(int SourceStateID, std::vector<int>* SuccIDV, std::vector<int>* CostV,
                          std::vector<EnvXXXLATAction_t*>* actionindV = NULL);

    int GetTrueCost(int parentID, int childID);
    bool isGoal(int id);

    /** \brief this function fill in Predecessor/Successor states of edges
     *         whose costs changed
     *         It takes in an array of cells whose traversability changed, and returns
     *         (in vector preds_of_changededgesIDV) the IDs of all states that have
     *         outgoing edges that go through the changed cells
     */
    void GetPredsofChangedEdges(std::vector<nav2dcell_t> const * changedcellsV,
                                        std::vector<int> *preds_of_changededgesIDV);

    /**
     * \brief same as GetPredsofChangedEdges, but returns successor states.
     *        Both functions need to be present for incremental search
     */
    void GetSuccsofChangedEdges(std::vector<nav2dcell_t> const * changedcellsV,
                                        std::vector<int> *succs_of_changededgesIDV);

    /**
     * \brief see comments on the same function in the parent class
     */
    void SetAllActionsandAllOutcomes(CMDPSTATE* state);

    /**
     * \brief see comments on the same function in the parent class
     */
    int GetFromToHeuristic(int FromStateID, int ToStateID);

    /**
     * \brief see comments on the same function in the parent class
     */
    int GetGoalHeuristic(int stateID);

    /**
     * \brief see comments on the same function in the parent class
     */
    int GetStartHeuristic(int stateID);

    /**
     * \brief see comments on the same function in the parent class
     */
    int SizeofCreatedEnv();

    /**
     * \brief see comments on the same function in the parent class
     */
    void PrintVars() { }

    const EnvXXXLATHashEntry_t* GetStateEntry(int state_id) const;

protected:
    //hash table of size x_size*y_size. Maps from coords to stateId
    int HashTableSize;
    std::vector<EnvXXXLATHashEntry_t*>* Coord2StateIDHashTable;
    //vector that maps from stateID to coords
    std::vector<EnvXXXLATHashEntry_t*> StateID2CoordTable;

    EnvXXXLATHashEntry_t** Coord2StateIDHashTable_lookup;

    unsigned int GETHASHBIN(unsigned int X, unsigned int Y, unsigned int Theta);

    EnvXXXLATHashEntry_t* GetHashEntry_hash(int X, int Y, int Theta);
    EnvXXXLATHashEntry_t* CreateNewHashEntry_hash(int X, int Y, int Theta);
    EnvXXXLATHashEntry_t* GetHashEntry_lookup(int X, int Y, int Theta);
    EnvXXXLATHashEntry_t* CreateNewHashEntry_lookup(int X, int Y, int Theta);

    //pointers to functions
    EnvXXXLATHashEntry_t* (EnvironmentXXXLAT::*GetHashEntry)(int X, int Y, int Theta);
    EnvXXXLATHashEntry_t* (EnvironmentXXXLAT::*CreateNewHashEntry)(int X, int Y, int Theta);

    void InitializeEnvironment();

    void PrintHashTableHist(FILE* fOut);
};

#endif
