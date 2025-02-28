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

#include <cmath>
#include <cstring>
#include <ctime>
#include <sbpl/discrete_space_information/environment_XXX.h>
#include <sbpl/utils/2Dgridsearch.h>
#include <sbpl/utils/key.h>
#include <sbpl/utils/mdp.h>
#include <sbpl/utils/mdpconfig.h>
// #include <fstream>

#if TIME_DEBUG
static clock_t time3_addallout = 0;
static clock_t time_gethash = 0;
static clock_t time_createhash = 0;
static clock_t time_getsuccs = 0;
#endif

static long int checks = 0;

#define XYTHETA2INDEX(X,Y,THETA) (THETA + X*EnvXXXCfg.NumThetaDirs + \
                                  Y*EnvXXXCfg.EnvWidth_c*EnvXXXCfg.NumThetaDirs)

EnvironmentXXXLATTICE::EnvironmentXXXLATTICE()
    : R0(0.56), // Distance from centre of robot to towbar
    F1(1.21), // Distance from towbar to front trailer wheel axle
    F2(0.87) // Distance from front trailer wheel axle to rear trailer wheel axle  
{
    EnvXXXCfg.obsthresh = EnvXXXLAT_DEFAULTOBSTHRESH;
    // the value that pretty much makes it disabled
    EnvXXXCfg.cost_inscribed_thresh = EnvXXXCfg.obsthresh;
    // the value that pretty much makes it disabled
    EnvXXXCfg.cost_possibly_circumscribed_thresh = -1;

    grid2Dsearchfromstart = NULL;
    grid2Dsearchfromgoal = NULL;
    bNeedtoRecomputeStartHeuristics = true;
    bNeedtoRecomputeGoalHeuristics = true;
    iteration = 0;
    bucketsize = 0; // fixed bucket size
    blocksize = 1;
    bUseNonUniformAngles = false;

    EnvXXXLAT.bInitialized = false;

    EnvXXXCfg.actionwidth = XXXLAT_DEFAULT_ACTIONWIDTH;

    EnvXXXCfg.NumThetaDirs = XXXLAT_THETADIRS;

    // no memory allocated in cfg yet
    EnvXXXCfg.Grid2D = NULL;
    EnvXXXCfg.ActionsV = NULL;
    EnvXXXCfg.PredActionsV = NULL;
}

EnvironmentXXXLATTICE::~EnvironmentXXXLATTICE()
{
    SBPL_PRINTF("destroying XYTHETALATTICE\n");
    if (grid2Dsearchfromstart != NULL) {
        delete grid2Dsearchfromstart;
    }
    grid2Dsearchfromstart = NULL;

    if (grid2Dsearchfromgoal != NULL) {
        delete grid2Dsearchfromgoal;
    }
    grid2Dsearchfromgoal = NULL;

    if (EnvXXXCfg.Grid2D != NULL) {
        for (int x = 0; x < EnvXXXCfg.EnvWidth_c; x++) {
            delete[] EnvXXXCfg.Grid2D[x];
        }
        delete[] EnvXXXCfg.Grid2D;
        EnvXXXCfg.Grid2D = NULL;
    }

    //delete actions
    if (EnvXXXCfg.ActionsV != NULL) {
        for (int tind = 0; tind < EnvXXXCfg.NumThetaDirs; tind++) {
            delete[] EnvXXXCfg.ActionsV[tind];
        }
        delete[] EnvXXXCfg.ActionsV;
        EnvXXXCfg.ActionsV = NULL;
    }
    if (EnvXXXCfg.PredActionsV != NULL) {
        delete[] EnvXXXCfg.PredActionsV;
        EnvXXXCfg.PredActionsV = NULL;
    }
}

static unsigned int inthash(unsigned int key)
{
    key += (key << 12);
    key ^= (key >> 22);
    key += (key << 4);
    key ^= (key >> 9);
    key += (key << 10);
    key ^= (key >> 2);
    key += (key << 7);
    key ^= (key >> 12);
    return key;
}

void EnvironmentXXXLATTICE::SetConfiguration(
    int width, int height, const unsigned char* mapdata,
    int startx, int starty, int starttheta, double starttheta1, double starttheta2,
    int goalx, int goaly, int goaltheta,
    double cellsize_m,
    double nominalvel_mpersecs,
    double timetoturn45degsinplace_secs,
    const std::vector<sbpl_2Dpt_t>& robot_perimeterV,
    const std::vector<sbpl_2Dpt_t>& trailer_perimeterV)
{
    EnvXXXCfg.EnvWidth_c = width;
    EnvXXXCfg.EnvHeight_c = height;
    EnvXXXCfg.StartX_c = startx;
    EnvXXXCfg.StartY_c = starty;
    EnvXXXCfg.StartTheta = starttheta;
    EnvXXXCfg.StartTheta1 = starttheta1;
    EnvXXXCfg.StartTheta2 = starttheta2;

    if (EnvXXXCfg.StartX_c < 0 ||
        EnvXXXCfg.StartX_c >= EnvXXXCfg.EnvWidth_c)
    {
        throw SBPL_Exception("ERROR: illegal start coordinates");
    }
    if (EnvXXXCfg.StartY_c < 0 ||
        EnvXXXCfg.StartY_c >= EnvXXXCfg.EnvHeight_c)
    {
        throw SBPL_Exception("ERROR: illegal start coordinates");
    }
    if (EnvXXXCfg.StartTheta < 0 ||
        EnvXXXCfg.StartTheta >= EnvXXXCfg.NumThetaDirs)
    {
        throw SBPL_Exception("ERROR: illegal start coordinates for theta");
    }

    EnvXXXCfg.EndX_c = goalx;
    EnvXXXCfg.EndY_c = goaly;
    EnvXXXCfg.EndTheta = goaltheta;

    if (EnvXXXCfg.EndX_c < 0 ||
        EnvXXXCfg.EndX_c >= EnvXXXCfg.EnvWidth_c)
    {
        throw SBPL_Exception("ERROR: illegal goal coordinates");
    }
    if (EnvXXXCfg.EndY_c < 0 ||
        EnvXXXCfg.EndY_c >= EnvXXXCfg.EnvHeight_c)
    {
        throw SBPL_Exception("ERROR: illegal goal coordinates");
    }
    if (EnvXXXCfg.EndTheta < 0 ||
        EnvXXXCfg.EndTheta >= EnvXXXCfg.NumThetaDirs)
    {
        throw SBPL_Exception("ERROR: illegal goal coordinates for theta");
    }

    EnvXXXCfg.FootprintPolygon = robot_perimeterV;
    EnvXXXCfg.TrailerPolygon = trailer_perimeterV;

    EnvXXXCfg.nominalvel_mpersecs = nominalvel_mpersecs;
    EnvXXXCfg.cellsize_m = cellsize_m;
    EnvXXXCfg.timetoturn45degsinplace_secs = timetoturn45degsinplace_secs;

    // unallocate the 2D environment
    if (EnvXXXCfg.Grid2D != NULL) {
        for (int x = 0; x < EnvXXXCfg.EnvWidth_c; x++) {
            delete[] EnvXXXCfg.Grid2D[x];
        }
        delete[] EnvXXXCfg.Grid2D;
        EnvXXXCfg.Grid2D = NULL;
    }

    // allocate the 2D environment
    EnvXXXCfg.Grid2D = new unsigned char*[EnvXXXCfg.EnvWidth_c];
    for (int x = 0; x < EnvXXXCfg.EnvWidth_c; x++) {
        EnvXXXCfg.Grid2D[x] = new unsigned char[EnvXXXCfg.EnvHeight_c];
    }

    // environment:
    if (0 == mapdata) {
        for (int y = 0; y < EnvXXXCfg.EnvHeight_c; y++) {
            for (int x = 0; x < EnvXXXCfg.EnvWidth_c; x++) {
                EnvXXXCfg.Grid2D[x][y] = 0;
            }
        }
    }
    else {
        for (int y = 0; y < EnvXXXCfg.EnvHeight_c; y++) {
            for (int x = 0; x < EnvXXXCfg.EnvWidth_c; x++) {
                EnvXXXCfg.Grid2D[x][y] = mapdata[x + y * width];
            }
        }
    }
}

void EnvironmentXXXLATTICE::ReadConfiguration(FILE* fCfg)
{
    // read in the configuration of environment and initialize
    // EnvXXXCfg structure
    char sTemp[1024], sTemp1[1024];
    int dTemp;
    int x, y;

    // discretization(cells)
    if (fscanf(fCfg, "%s", sTemp) != 1) {
        throw SBPL_Exception("ERROR: ran out of env file early (discretization)");
    }
    strcpy(sTemp1, "discretization(cells):");
    if (strcmp(sTemp1, sTemp) != 0) {
        std::stringstream ss;
        ss << "ERROR: configuration file has incorrect format (discretization)" <<
                " Expected " << sTemp1 << " got " << sTemp;
        throw SBPL_Exception(ss.str());
    }
    if (fscanf(fCfg, "%s", sTemp) != 1) {
        throw SBPL_Exception("ERROR: ran out of env file early (discretization)");
    }
    EnvXXXCfg.EnvWidth_c = atoi(sTemp);
    if (fscanf(fCfg, "%s", sTemp) != 1) {
        throw SBPL_Exception("ERROR: ran out of env file early (discretization)");
    }
    EnvXXXCfg.EnvHeight_c = atoi(sTemp);

    // Scan for optional NumThetaDirs parameter. Check for following obsthresh.
    if (fscanf(fCfg, "%s", sTemp) != 1) {
        throw SBPL_Exception("ERROR: ran out of env file early");
    }
    strcpy(sTemp1, "NumThetaDirs:");
    if (strcmp(sTemp1, sTemp) != 0) {
        // optional NumThetaDirs not available; default is XXXLAT_THETADIRS (16)
        strcpy(sTemp1, "obsthresh:");
        if (strcmp(sTemp1, sTemp) != 0) {
            std::stringstream ss;
            ss << "ERROR: configuration file has incorrect format" <<
                    " Expected " << sTemp1 << " got " << sTemp;
            throw SBPL_Exception(ss.str());
        }
        else {
            EnvXXXCfg.NumThetaDirs = XXXLAT_THETADIRS;
        }
    }
    else {
        if (fscanf(fCfg, "%s", sTemp) != 1) {
            throw SBPL_Exception("ERROR: ran out of env file early (NumThetaDirs)");
        }
        EnvXXXCfg.NumThetaDirs = atoi(sTemp);

        //obsthresh:
        if (fscanf(fCfg, "%s", sTemp) != 1) {
            throw SBPL_Exception("ERROR: ran out of env file early (obsthresh)");
        }
        strcpy(sTemp1, "obsthresh:");
        if (strcmp(sTemp1, sTemp) != 0) {
            std::stringstream ss;
            ss << "ERROR: configuration file has incorrect format" <<
                " Expected " << sTemp1 << " got " << sTemp <<
                " see existing examples of env files for the right format of heading";
            throw SBPL_Exception(ss.str());
        }
    }

    // obsthresh
    if (fscanf(fCfg, "%s", sTemp) != 1) {
        throw SBPL_Exception("ERROR: ran out of env file early");
    }
    EnvXXXCfg.obsthresh = atoi(sTemp);
    SBPL_PRINTF("obsthresh = %d\n", EnvXXXCfg.obsthresh);

    //cost_inscribed_thresh:
    if (fscanf(fCfg, "%s", sTemp) != 1) {
        throw SBPL_Exception("ERROR: ran out of env file early");
    }
    strcpy(sTemp1, "cost_inscribed_thresh:");
    if (strcmp(sTemp1, sTemp) != 0) {
        std::stringstream ss;
        ss << "ERROR: configuration file has incorrect format" <<
                " Expected " << sTemp1 << " got " << sTemp <<
                " see existing examples of env files for the right format of heading";
        throw SBPL_Exception(ss.str());
    }
    if (fscanf(fCfg, "%s", sTemp) != 1) {
        throw SBPL_Exception("ERROR: ran out of env file early");
    }
    EnvXXXCfg.cost_inscribed_thresh = atoi(sTemp);
    SBPL_PRINTF("cost_inscribed_thresh = %d\n", EnvXXXCfg.cost_inscribed_thresh);

    //cost_possibly_circumscribed_thresh:
    if (fscanf(fCfg, "%s", sTemp) != 1) {
        throw SBPL_Exception("ERROR: ran out of env file early");
    }
    strcpy(sTemp1, "cost_possibly_circumscribed_thresh:");
    if (strcmp(sTemp1, sTemp) != 0) {
        std::stringstream ss;
        ss << "ERROR: configuration file has incorrect format" <<
                " Expected " << sTemp1 << " got " << sTemp <<
                " see existing examples of env files for the right format of heading";
        throw SBPL_Exception(ss.str());
    }
    if (fscanf(fCfg, "%s", sTemp) != 1) {
        throw SBPL_Exception("ERROR: ran out of env file early");
    }
    EnvXXXCfg.cost_possibly_circumscribed_thresh = atoi(sTemp);
    SBPL_PRINTF("cost_possibly_circumscribed_thresh = %d\n", EnvXXXCfg.cost_possibly_circumscribed_thresh);

    //cellsize
    if (fscanf(fCfg, "%s", sTemp) != 1) {
        throw SBPL_Exception("ERROR: ran out of env file early");
    }
    strcpy(sTemp1, "cellsize(meters):");
    if (strcmp(sTemp1, sTemp) != 0) {
        std::stringstream ss;
        ss << "ERROR: configuration file has incorrect format" <<
                " Expected " << sTemp1 << " got " << sTemp;
        throw SBPL_Exception(ss.str());
    }
    if (fscanf(fCfg, "%s", sTemp) != 1) {
        throw SBPL_Exception("ERROR: ran out of env file early");
    }
    EnvXXXCfg.cellsize_m = atof(sTemp);

    //speeds
    if (fscanf(fCfg, "%s", sTemp) != 1) {
        throw SBPL_Exception("ERROR: ran out of env file early");
    }
    strcpy(sTemp1, "nominalvel(mpersecs):");
    if (strcmp(sTemp1, sTemp) != 0) {
        std::stringstream ss;
        ss << "ERROR: configuration file has incorrect format" <<
            " Expected " << sTemp1 << " got " << sTemp;
        throw SBPL_Exception(ss.str());
    }
    if (fscanf(fCfg, "%s", sTemp) != 1) {
        throw SBPL_Exception("ERROR: ran out of env file early");
    }
    EnvXXXCfg.nominalvel_mpersecs = atof(sTemp);
    if (fscanf(fCfg, "%s", sTemp) != 1) {
        throw SBPL_Exception("ERROR: ran out of env file early");
    }
    strcpy(sTemp1, "timetoturn45degsinplace(secs):");
    if (strcmp(sTemp1, sTemp) != 0) {
        std::stringstream ss;
        ss << "ERROR: configuration file has incorrect format" <<
                " Expected " << sTemp1 << " got " << sTemp;
        throw SBPL_Exception(ss.str());
    }
    if (fscanf(fCfg, "%s", sTemp) != 1) {
        throw SBPL_Exception("ERROR: ran out of env file early");
    }
    EnvXXXCfg.timetoturn45degsinplace_secs = atof(sTemp);

    // start(meters,rads):
    if (fscanf(fCfg, "%s", sTemp) != 1) {
        throw SBPL_Exception("ERROR: ran out of env file early");
    }
    if (fscanf(fCfg, "%s", sTemp) != 1) {
        throw SBPL_Exception("ERROR: ran out of env file early");
    }
    EnvXXXCfg.StartX_c = CONTXY2DISC(atof(sTemp), EnvXXXCfg.cellsize_m);
    if (fscanf(fCfg, "%s", sTemp) != 1) {
        throw SBPL_Exception("ERROR: ran out of env file early");
    }
    EnvXXXCfg.StartY_c = CONTXY2DISC(atof(sTemp), EnvXXXCfg.cellsize_m);
    if (fscanf(fCfg, "%s", sTemp) != 1) {
        throw SBPL_Exception("ERROR: ran out of env file early");
    }

    EnvXXXCfg.StartTheta_rad = atof(sTemp);

    if (EnvXXXCfg.StartX_c < 0 ||
        EnvXXXCfg.StartX_c >= EnvXXXCfg.EnvWidth_c)
    {
        throw SBPL_Exception("ERROR: illegal start coordinates");
    }
    if (EnvXXXCfg.StartY_c < 0 ||
        EnvXXXCfg.StartY_c >= EnvXXXCfg.EnvHeight_c)
    {
        throw SBPL_Exception("ERROR: illegal start coordinates");
    }

    // end(meters,rads):
    if (fscanf(fCfg, "%s", sTemp) != 1) {
        throw SBPL_Exception("ERROR: ran out of env file early");
    }
    if (fscanf(fCfg, "%s", sTemp) != 1) {
        throw SBPL_Exception("ERROR: ran out of env file early");
    }
    EnvXXXCfg.EndX_c = CONTXY2DISC(atof(sTemp), EnvXXXCfg.cellsize_m);
    if (fscanf(fCfg, "%s", sTemp) != 1) {
        throw SBPL_Exception("ERROR: ran out of env file early");
    }
    EnvXXXCfg.EndY_c = CONTXY2DISC(atof(sTemp), EnvXXXCfg.cellsize_m);
    if (fscanf(fCfg, "%s", sTemp) != 1) {
        throw SBPL_Exception("ERROR: ran out of env file early");
    }

    EnvXXXCfg.EndTheta_rad = atof(sTemp);

    if (EnvXXXCfg.EndX_c < 0 ||
        EnvXXXCfg.EndX_c >= EnvXXXCfg.EnvWidth_c)
    {
        throw SBPL_Exception("ERROR: illegal end coordinates");
    }
    if (EnvXXXCfg.EndY_c < 0 ||
        EnvXXXCfg.EndY_c >= EnvXXXCfg.EnvHeight_c)
    {
        throw SBPL_Exception("ERROR: illegal end coordinates");
    }

    // unallocate the 2d environment
    if (EnvXXXCfg.Grid2D != NULL) {
        for (x = 0; x < EnvXXXCfg.EnvWidth_c; x++) {
            delete[] EnvXXXCfg.Grid2D[x];
        }
        delete[] EnvXXXCfg.Grid2D;
        EnvXXXCfg.Grid2D = NULL;
    }

    // allocate the 2D environment
    EnvXXXCfg.Grid2D = new unsigned char*[EnvXXXCfg.EnvWidth_c];
    for (x = 0; x < EnvXXXCfg.EnvWidth_c; x++) {
        EnvXXXCfg.Grid2D[x] = new unsigned char[EnvXXXCfg.EnvHeight_c];
    }

    // environment:
    if (fscanf(fCfg, "%s", sTemp) != 1) {
        throw SBPL_Exception("ERROR: ran out of env file early");
    }
    for (y = 0; y < EnvXXXCfg.EnvHeight_c; y++) {
        for (x = 0; x < EnvXXXCfg.EnvWidth_c; x++) {
            if (fscanf(fCfg, "%d", &dTemp) != 1) {
                throw SBPL_Exception("ERROR: incorrect format of config file");
            }
            EnvXXXCfg.Grid2D[x][y] = dTemp;
        }
    }
}

bool EnvironmentXXXLATTICE::ReadinCell(
    sbpl_xy_theta_cell_t* cell,
    FILE* fIn)
{
    char sTemp[60];

    if (fscanf(fIn, "%s", sTemp) == 0) {
        return false;
    }
    cell->x = atoi(sTemp);
    if (fscanf(fIn, "%s", sTemp) == 0) {
        return false;
    }
    cell->y = atoi(sTemp);
    if (fscanf(fIn, "%s", sTemp) == 0) {
        return false;
    }
    cell->theta = atoi(sTemp);

    // normalize the angle
    cell->theta = normalizeDiscAngle(cell->theta);
    // cell->theta = NORMALIZEDISCTHETA(cell->theta, EnvXXXCfg.NumThetaDirs);

    return true;
}

bool EnvironmentXXXLATTICE::ReadinPose(
    sbpl_xy_theta_pt_t* pose,
    FILE* fIn)
{
    char sTemp[60];

    if (fscanf(fIn, "%s", sTemp) == 0) {
        return false;
    }
    pose->x = atof(sTemp);
    if (fscanf(fIn, "%s", sTemp) == 0) {
        return false;
    }
    pose->y = atof(sTemp);
    if (fscanf(fIn, "%s", sTemp) == 0) {
        return false;
    }
    pose->theta = atof(sTemp);

    pose->theta = normalizeAngle(pose->theta);

    return true;
}

int EnvironmentXXXLATTICE::normalizeDiscAngle(int theta) const
{
    if (bUseNonUniformAngles) {
        if (theta < 0) {
            theta += EnvXXXCfg.NumThetaDirs;
        }
        if (theta >= EnvXXXCfg.NumThetaDirs) {
            theta -= EnvXXXCfg.NumThetaDirs;
        }
    }
    else {
        theta = NORMALIZEDISCTHETA(theta, EnvXXXCfg.NumThetaDirs);
    }
    return theta;
}

double EnvironmentXXXLATTICE::DiscTheta2ContNew(int theta) const
{
    if (bUseNonUniformAngles) {
        return DiscTheta2ContFromSet(theta);
    }
    else {
        return DiscTheta2Cont(theta, EnvXXXCfg.NumThetaDirs);
    }
}

int EnvironmentXXXLATTICE::ContTheta2DiscNew(double theta) const
{
    if (bUseNonUniformAngles) {
        return ContTheta2DiscFromSet(theta);
    }
    else {
        return ContTheta2Disc(theta, EnvXXXCfg.NumThetaDirs);
    }
}

double EnvironmentXXXLATTICE::DiscTheta2ContFromSet(int theta) const
{
    theta = normalizeDiscAngle(theta);

    // ThetaDirs should contain extra angle (2PI) for overlap
    if (EnvXXXCfg.NumThetaDirs >= (int)EnvXXXCfg.ThetaDirs.size()) {
        throw SBPL_Exception("ERROR: list of bin angles are not properly set to use function DiscTheta2ConfFromSet");
    }

    if (theta > EnvXXXCfg.NumThetaDirs || theta < 0) {
        std::stringstream ss;
        ss << "ERROR: discrete value theta " << theta << " out of range";
        throw SBPL_Exception(ss.str());
    }
    return EnvXXXCfg.ThetaDirs[theta];
}

int EnvironmentXXXLATTICE::ContTheta2DiscFromSet(double theta) const
{
    theta = normalizeAngle(theta);
    // ThetaDirs should contain extra angle (2PI) for overlap
    if (EnvXXXCfg.NumThetaDirs >= (int) EnvXXXCfg.ThetaDirs.size()) {
        throw SBPL_Exception("ERROR: list of bin angles are not properly set to use function ContTheta2DiscFromSet");
    }

    int lower_bound_ind = -1;
    int upper_bound_ind = -1;
    for (int i = 1; i < (int) EnvXXXCfg.ThetaDirs.size(); i++) {
        if ((EnvXXXCfg.ThetaDirs[i]) >= theta) {
            lower_bound_ind = i - 1;
            upper_bound_ind = i;
            break;
        }
    }

    // Critical error if could not find bin location from given angle
    if (lower_bound_ind == -1) {
        std::stringstream ss;
        ss << "ERROR: unable to find bin index for angle " << theta;
        throw SBPL_Exception(ss.str());
    }

    // Get closest angle of two
    double angle_low = EnvXXXCfg.ThetaDirs[lower_bound_ind];
    double angle_up = EnvXXXCfg.ThetaDirs[upper_bound_ind];
    double diff_low = fabs(theta - angle_low);
    double diff_up = fabs(theta - angle_up);

    if (diff_low < diff_up) {
        return lower_bound_ind;
    }
    else {
        // Wrap upper bound index around when it reaches last index (assumed to be 2PI)
        if (upper_bound_ind == EnvXXXCfg.NumThetaDirs) {
            upper_bound_ind = 0;
        }
        return upper_bound_ind;
    }
}

bool EnvironmentXXXLATTICE::ReadinMotionPrimitive(
    SBPL_xxx_mprimitive* pMotPrim,
    FILE* fIn)
{
    char sTemp[1024];
    int dTemp;
    char sExpected[1024];
    int numofIntermPoses;
    float fTemp;

    // read in actionID
    strcpy(sExpected, "primID:");
    if (fscanf(fIn, "%s", sTemp) == 0) {
        return false;
    }
    if (strcmp(sTemp, sExpected) != 0) {
        SBPL_ERROR("ERROR: expected %s but got %s\n", sExpected, sTemp);
        fflush(stdout);
        return false;
    }
    if (fscanf(fIn, "%d", &pMotPrim->motprimID) != 1) {
        return false;
    }

    // read in start angle
    strcpy(sExpected, "startangle_c:");
    if (fscanf(fIn, "%s", sTemp) == 0) {
        return false;
    }
    if (strcmp(sTemp, sExpected) != 0) {
        SBPL_ERROR("ERROR: expected %s but got %s\n", sExpected, sTemp);
        return false;
    }
    if (fscanf(fIn, "%d", &dTemp) == 0) {
        SBPL_ERROR("ERROR reading startangle\n");
        return false;
    }
    pMotPrim->starttheta_c = dTemp;

    // read in end pose
    strcpy(sExpected, "endpose_c:");
    if (fscanf(fIn, "%s", sTemp) == 0) {
        return false;
    }
    if (strcmp(sTemp, sExpected) != 0) {
        SBPL_ERROR("ERROR: expected %s but got %s\n", sExpected, sTemp);
        return false;
    }

    if (ReadinCell(&pMotPrim->endcell, fIn) == false) {
        SBPL_ERROR("ERROR: failed to read in endsearchpose\n");
        return false;
    }

    // read in action cost
    strcpy(sExpected, "additionalactioncostmult:");
    if (fscanf(fIn, "%s", sTemp) == 0) {
        return false;
    }
    if (strcmp(sTemp, sExpected) != 0) {
        SBPL_ERROR("ERROR: expected %s but got %s\n", sExpected, sTemp);
        return false;
    }
    if (fscanf(fIn, "%d", &dTemp) != 1) {
        return false;
    }
    pMotPrim->additionalactioncostmult = dTemp;

    if (bUseNonUniformAngles) {
        // read in action turning radius
        strcpy(sExpected, "turning_radius:");
        if (fscanf(fIn, "%s", sTemp) == 0) {
            return false;
        }
        if (strcmp(sTemp, sExpected) != 0) {
            SBPL_ERROR("ERROR: expected %s but got %s\n", sExpected, sTemp);
            return false;
        }
        if (fscanf(fIn, "%f", &fTemp) != 1) {
            return false;
        }
        pMotPrim->turning_radius = fTemp;
    }

    // read in intermediate poses
    strcpy(sExpected, "intermediateposes:");
    if (fscanf(fIn, "%s", sTemp) == 0) {
        return false;
    }
    if (strcmp(sTemp, sExpected) != 0) {
        SBPL_ERROR("ERROR: expected %s but got %s\n", sExpected, sTemp);
        return false;
    }
    if (fscanf(fIn, "%d", &numofIntermPoses) != 1) {
        return false;
    }
    // all intermposes should be with respect to 0,0 as starting pose since it
    // will be added later and should be done after the action is rotated by
    // initial orientation
    for (int i = 0; i < numofIntermPoses; i++) {
        sbpl_xy_theta_pt_t intermpose;
        if (ReadinPose(&intermpose, fIn) == false) {
            SBPL_ERROR("ERROR: failed to read in intermediate poses\n");
            return false;
        }
        pMotPrim->intermptV.push_back(intermpose);
    }

    // Check that the last pose of the motion matches (within lattice
    // resolution) the designated end pose of the primitive
    sbpl_xy_theta_pt_t sourcepose;
    sourcepose.x = DISCXY2CONT(0, EnvXXXCfg.cellsize_m);
    sourcepose.y = DISCXY2CONT(0, EnvXXXCfg.cellsize_m);
    sourcepose.theta = DiscTheta2ContNew(pMotPrim->starttheta_c);
    double mp_endx_m = sourcepose.x + pMotPrim->intermptV[pMotPrim->intermptV.size() - 1].x;
    double mp_endy_m = sourcepose.y + pMotPrim->intermptV[pMotPrim->intermptV.size() - 1].y;
    double mp_endtheta_rad = pMotPrim->intermptV[pMotPrim->intermptV.size() - 1].theta;

    int endtheta_c;
    int endx_c = CONTXY2DISC(mp_endx_m, EnvXXXCfg.cellsize_m);
    int endy_c = CONTXY2DISC(mp_endy_m, EnvXXXCfg.cellsize_m);
    endtheta_c = ContTheta2DiscNew(mp_endtheta_rad);
    if (endx_c != pMotPrim->endcell.x ||
        endy_c != pMotPrim->endcell.y ||
        endtheta_c != pMotPrim->endcell.theta)
    {
        SBPL_ERROR( "ERROR: incorrect primitive %d with startangle=%d "
                   "last interm point %f %f %f does not match end pose %d %d %d\n",
                   pMotPrim->motprimID, pMotPrim->starttheta_c,
                   pMotPrim->intermptV[pMotPrim->intermptV.size() - 1].x,
                   pMotPrim->intermptV[pMotPrim->intermptV.size() - 1].y,
                   pMotPrim->intermptV[pMotPrim->intermptV.size() - 1].theta,
                   pMotPrim->endcell.x, pMotPrim->endcell.y,
                   pMotPrim->endcell.theta);
        SBPL_FFLUSH(stdout);
        return false;
    }

    return true;
}

bool EnvironmentXXXLATTICE::ReadMotionPrimitives(FILE* fMotPrims)
{
    char sTemp[1024], sExpected[1024];
    float fTemp;
    int dTemp;
    int totalNumofActions = 0;

    SBPL_INFO("Reading in motion primitives...");
    fflush(stdout);

    //read in the resolution
    strcpy(sExpected, "resolution_m:");
    if (fscanf(fMotPrims, "%s", sTemp) == 0) {
        return false;
    }
    if (strcmp(sTemp, sExpected) != 0) {
        SBPL_ERROR("ERROR: expected %s but got %s\n", sExpected, sTemp);
        fflush(stdout);
        return false;
    }
    if (fscanf(fMotPrims, "%f", &fTemp) == 0) {
        return false;
    }
    if (fabs(fTemp - EnvXXXCfg.cellsize_m) > ERR_EPS) {
        SBPL_ERROR("ERROR: invalid resolution %f (instead of %f) in the dynamics file\n", fTemp, EnvXXXCfg.cellsize_m);
        fflush(stdout);
        return false;
    }
    SBPL_INFO("resolution_m: %f\n", fTemp);

    if (fscanf(fMotPrims, "%s", sTemp) == 0) {
        return false;
    }
    SBPL_INFO("sTemp: %s\n", sTemp);
    if (strncmp(sTemp, "min_turning_radius_m:", 21) == 0) {
        bUseNonUniformAngles = true;
    }
    SBPL_INFO("bUseNonUniformAngles = %d", bUseNonUniformAngles);

    if (bUseNonUniformAngles) {
        float min_turn_rad;
        strcpy(sExpected, "min_turning_radius_m:");
        if (strcmp(sTemp, sExpected) != 0) {
            SBPL_ERROR("ERROR: expected %s but got %s\n", sExpected, sTemp);
            fflush(stdout);
            return false;
        }
        if (fscanf(fMotPrims, "%f", &min_turn_rad) == 0) {
            return false;
        }
        SBPL_PRINTF("min_turn_rad: %f\n", min_turn_rad);
        fflush(stdout);
        if (fscanf(fMotPrims, "%s", sTemp) == 0) {
            return false;
        }
    }

    // read in the angular resolution
    strcpy(sExpected, "numberofangles:");
    if (strcmp(sTemp, sExpected) != 0) {
       SBPL_ERROR("ERROR: expected %s but got %s\n", sExpected, sTemp);
       return false;
    }
    if (fscanf(fMotPrims, "%d", &dTemp) == 0) {
        return false;
    }
    if (dTemp != EnvXXXCfg.NumThetaDirs) {
        SBPL_ERROR("ERROR: invalid angular resolution %d angles (instead of %d angles) in the motion primitives file\n", dTemp, EnvXXXCfg.NumThetaDirs);
        return false;
    }
    SBPL_PRINTF("numberofangles: %d\n", dTemp);
    EnvXXXCfg.NumThetaDirs = dTemp;

    if (bUseNonUniformAngles) {
        // read in angles
        EnvXXXCfg.ThetaDirs.clear();
        for (int i = 0; i < EnvXXXCfg.NumThetaDirs; i++)
        {
            std::ostringstream string_angle_index;
            string_angle_index << i;
            std::string angle_string = "angle:" + string_angle_index.str();

            float angle;
            strcpy(sExpected, angle_string.c_str());
            if (fscanf(fMotPrims, "%s", sTemp) == 0) {
                return false;
            }
            if (strcmp(sTemp, sExpected) != 0) {
                SBPL_ERROR("ERROR: expected %s but got %s\n", sExpected, sTemp);
                return false;
            }
            if (fscanf(fMotPrims, "%f", &angle) == 0) {
                return false;
            }
            SBPL_PRINTF("%s %f\n", angle_string.c_str(), angle);
            EnvXXXCfg.ThetaDirs.push_back(angle);
        }
        EnvXXXCfg.ThetaDirs.push_back(2.0 * M_PI); // Add 2 PI at end for overlap
    }

    // read in the total number of actions
    strcpy(sExpected, "totalnumberofprimitives:");
    if (fscanf(fMotPrims, "%s", sTemp) == 0) {
        return false;
    }
    if (strcmp(sTemp, sExpected) != 0) {
        SBPL_ERROR("ERROR: expected %s but got %s\n", sExpected, sTemp);
        return false;
    }
    if (fscanf(fMotPrims, "%d", &totalNumofActions) == 0) {
        return false;
    }
    SBPL_PRINTF("totalnumberofprimitives: %d\n", totalNumofActions);

    // Read in motion primitive for each action
    for (int i = 0; i < totalNumofActions; i++) {
        SBPL_xxx_mprimitive motprim;

        if (!EnvironmentXXXLATTICE::ReadinMotionPrimitive(&motprim, fMotPrims)) {
            return false;
        }

        EnvXXXCfg.mprimV.push_back(motprim);
    }
    SBPL_PRINTF("done");
    SBPL_FFLUSH(stdout);
    return true;
}

void EnvironmentXXXLATTICE::ComputeReplanningDataforAction(
    EnvXXXLATAction_t* action)
{
    int j;

    // iterate over all the cells involved in the action
    sbpl_xy_theta_cell_t startcell3d, endcell3d;
    for (int i = 0; i < (int)action->intersectingcellsV.size(); i++) {
        // compute the translated affected search Pose - what state has an
        // outgoing action whose intersecting cell is at 0,0
        startcell3d.theta = action->starttheta;
        startcell3d.x = -action->intersectingcellsV.at(i).x;
        startcell3d.y = -action->intersectingcellsV.at(i).y;

        // compute the translated affected search Pose - what state has an
        // incoming action whose intersecting cell is at 0,0
        if (bUseNonUniformAngles) {
            endcell3d.theta = normalizeDiscAngle(action->endtheta);
        }
        else {
            endcell3d.theta = NORMALIZEDISCTHETA(action->endtheta, EnvXXXCfg.NumThetaDirs);
        }
        endcell3d.x = startcell3d.x + action->dX;
        endcell3d.y = startcell3d.y + action->dY;

        //store the cells if not already there
        for (j = 0; j < (int)affectedsuccstatesV.size(); j++) {
            if (affectedsuccstatesV.at(j) == endcell3d) {
                break;
            }
        }
        if (j == (int)affectedsuccstatesV.size()) {
            affectedsuccstatesV.push_back(endcell3d);
        }

        for (j = 0; j < (int)affectedpredstatesV.size(); j++) {
            if (affectedpredstatesV.at(j) == startcell3d) {
                break;
            }
        }
        if (j == (int)affectedpredstatesV.size()) {
            affectedpredstatesV.push_back(startcell3d);
        }
    } // over intersecting cells

    // add the centers since with h2d we are using these in cost computations
    // ---intersecting cell = origin
    // compute the translated affected search Pose - what state has an outgoing
    // action whose intersecting cell is at 0,0
    startcell3d.theta = action->starttheta;
    startcell3d.x = -0;
    startcell3d.y = -0;

    // compute the translated affected search Pose - what state has an incoming
    // action whose intersecting cell is at 0,0
    if (bUseNonUniformAngles) {
        endcell3d.theta = normalizeDiscAngle(action->endtheta);
    }
    else {
        endcell3d.theta = NORMALIZEDISCTHETA(action->endtheta, EnvXXXCfg.NumThetaDirs);
    }
    endcell3d.x = startcell3d.x + action->dX;
    endcell3d.y = startcell3d.y + action->dY;

    //store the cells if not already there
    for (j = 0; j < (int)affectedsuccstatesV.size(); j++) {
        if (affectedsuccstatesV.at(j) == endcell3d) {
            break;
        }
    }
    if (j == (int)affectedsuccstatesV.size()) {
        affectedsuccstatesV.push_back(endcell3d);
    }

    for (j = 0; j < (int)affectedpredstatesV.size(); j++) {
        if (affectedpredstatesV.at(j) == startcell3d) {
            break;
        }
    }
    if (j == (int)affectedpredstatesV.size()) {
        affectedpredstatesV.push_back(startcell3d);
    }

    //---intersecting cell = outcome state
    // compute the translated affected search Pose - what state has an outgoing
    // action whose intersecting cell is at 0,0
    startcell3d.theta = action->starttheta;
    startcell3d.x = -action->dX;
    startcell3d.y = -action->dY;

    // compute the translated affected search Pose - what state has an incoming
    // action whose intersecting cell is at 0,0
    if (bUseNonUniformAngles) {
        endcell3d.theta = normalizeDiscAngle(action->endtheta);
    }
    else {
        endcell3d.theta = NORMALIZEDISCTHETA(action->endtheta, EnvXXXCfg.NumThetaDirs);
    }
    endcell3d.x = startcell3d.x + action->dX;
    endcell3d.y = startcell3d.y + action->dY;

    for (j = 0; j < (int)affectedsuccstatesV.size(); j++) {
        if (affectedsuccstatesV.at(j) == endcell3d) {
            break;
        }
    }
    if (j == (int)affectedsuccstatesV.size()) {
        affectedsuccstatesV.push_back(endcell3d);
    }

    for (j = 0; j < (int)affectedpredstatesV.size(); j++) {
        if (affectedpredstatesV.at(j) == startcell3d) {
            break;
        }
    }
    if (j == (int)affectedpredstatesV.size()) {
        affectedpredstatesV.push_back(startcell3d);
    }
}

// computes all the 3D states whose outgoing actions are potentially affected
// when cell (0,0) changes its status it also does the same for the 3D states
// whose incoming actions are potentially affected when cell (0,0) changes its
// status
void EnvironmentXXXLATTICE::ComputeReplanningData()
{
    // iterate over all actions
    // orientations
    for (int tind = 0; tind < EnvXXXCfg.NumThetaDirs; tind++) {
        // actions
        for (int aind = 0; aind < EnvXXXCfg.actionwidth; aind++) {
            // compute replanning data for this action
            ComputeReplanningDataforAction(&EnvXXXCfg.ActionsV[tind][aind]);
        }
    }
}

// here motionprimitivevector contains actions only for 0 angle
void EnvironmentXXXLATTICE::PrecomputeActionswithBaseMotionPrimitive(
    std::vector<SBPL_xxx_mprimitive>* motionprimitiveV)
{
    SBPL_PRINTF("Pre-computing action data using base motion primitives...\n");
    EnvXXXCfg.ActionsV = new EnvXXXLATAction_t*[EnvXXXCfg.NumThetaDirs];
    EnvXXXCfg.PredActionsV = new std::vector<EnvXXXLATAction_t*> [EnvXXXCfg.NumThetaDirs];
    std::vector<sbpl_2Dcell_t> footprint;

    //iterate over source angles
    for (int tind = 0; tind < EnvXXXCfg.NumThetaDirs; tind++) {
        SBPL_PRINTF("pre-computing for angle %d out of %d angles\n", tind, EnvXXXCfg.NumThetaDirs);
        EnvXXXCfg.ActionsV[tind] = new EnvXXXLATAction_t[motionprimitiveV->size()];

        //compute sourcepose
        sbpl_xy_theta_pt_t sourcepose;
        sourcepose.x = DISCXY2CONT(0, EnvXXXCfg.cellsize_m);
        sourcepose.y = DISCXY2CONT(0, EnvXXXCfg.cellsize_m);
    sourcepose.theta = DiscTheta2ContNew(tind);
    //sourcepose.theta = DiscTheta2Cont(tind, EnvXXXCfg.NumThetaDirs);

        //iterate over motion primitives
        for (size_t aind = 0; aind < motionprimitiveV->size(); aind++) {
            EnvXXXCfg.ActionsV[tind][aind].aind = aind;
            EnvXXXCfg.ActionsV[tind][aind].starttheta = tind;
            double mp_endx_m = motionprimitiveV->at(aind).intermptV[motionprimitiveV->at(aind).intermptV.size() - 1].x;
            double mp_endy_m = motionprimitiveV->at(aind).intermptV[motionprimitiveV->at(aind).intermptV.size() - 1].y;
            double mp_endtheta_rad =
                    motionprimitiveV->at(aind).intermptV[motionprimitiveV->at(aind).intermptV.size() - 1].theta;

            double endx = sourcepose.x + (mp_endx_m * cos(sourcepose.theta) - mp_endy_m * sin(sourcepose.theta));
            double endy = sourcepose.y + (mp_endx_m * sin(sourcepose.theta) + mp_endy_m * cos(sourcepose.theta));

            int endx_c = CONTXY2DISC(endx, EnvXXXCfg.cellsize_m);
            int endy_c = CONTXY2DISC(endy, EnvXXXCfg.cellsize_m);
        EnvXXXCfg.ActionsV[tind][aind].endtheta = ContTheta2DiscNew(mp_endtheta_rad + sourcepose.theta);
        //EnvXXXCfg.ActionsV[tind][aind].endtheta = ContTheta2Disc(mp_endtheta_rad + sourcepose.theta, EnvXXXCfg.NumThetaDirs);
            EnvXXXCfg.ActionsV[tind][aind].dX = endx_c;
            EnvXXXCfg.ActionsV[tind][aind].dY = endy_c;

            if (EnvXXXCfg.ActionsV[tind][aind].dY != 0 || EnvXXXCfg.ActionsV[tind][aind].dX != 0)
                EnvXXXCfg.ActionsV[tind][aind].cost = (int)(ceil(XXXLAT_COSTMULT_MTOMM
                    * EnvXXXCfg.cellsize_m / EnvXXXCfg.nominalvel_mpersecs
                    * sqrt((double)(EnvXXXCfg.ActionsV[tind][aind].dX
                        * EnvXXXCfg.ActionsV[tind][aind].dX + EnvXXXCfg.ActionsV[tind][aind].dY
                        * EnvXXXCfg.ActionsV[tind][aind].dY))));
            else
                //cost of turn in place
                EnvXXXCfg.ActionsV[tind][aind].cost = (int)(XXXLAT_COSTMULT_MTOMM
                    * EnvXXXCfg.timetoturn45degsinplace_secs
                    * fabs(computeMinUnsignedAngleDiff(mp_endtheta_rad, 0)) / (PI_CONST / 4.0));

            //compute and store interm points as well as intersecting cells
            EnvXXXCfg.ActionsV[tind][aind].intersectingcellsV.clear();
            EnvXXXCfg.ActionsV[tind][aind].intermptV.clear();
            EnvXXXCfg.ActionsV[tind][aind].interm3DcellsV.clear();
            sbpl_xy_theta_cell_t previnterm3Dcell;
            previnterm3Dcell.theta = previnterm3Dcell.x = previnterm3Dcell.y = 0;

            for (int pind = 0; pind < (int)motionprimitiveV->at(aind).intermptV.size(); pind++) {
                sbpl_xy_theta_pt_t intermpt = motionprimitiveV->at(aind).intermptV[pind];

                //rotate it appropriately
                double rotx = intermpt.x * cos(sourcepose.theta) - intermpt.y * sin(sourcepose.theta);
                double roty = intermpt.x * sin(sourcepose.theta) + intermpt.y * cos(sourcepose.theta);
                intermpt.x = rotx;
                intermpt.y = roty;
                intermpt.theta = normalizeAngle(sourcepose.theta + intermpt.theta);

                //store it (they are with reference to 0,0,stattheta (not
                //sourcepose.x,sourcepose.y,starttheta (that is, half-bin))
                EnvXXXCfg.ActionsV[tind][aind].intermptV.push_back(intermpt);
            }
            //now compute the intersecting cells for this motion (including ignoring the source footprint)
            get_2d_motion_cells(EnvXXXCfg.FootprintPolygon,
                                EnvXXXCfg.ActionsV[tind][aind].intermptV,
                                &EnvXXXCfg.ActionsV[tind][aind].intersectingcellsV,
                                EnvXXXCfg.cellsize_m);

#if DEBUG
            SBPL_FPRINTF(fDeb,
                         "action tind=%d aind=%d: dX=%d dY=%d endtheta=%d (%.2f degs -> %.2f degs) "
                         "cost=%d (mprim: %.2f %.2f %.2f)\n",
                         tind,
                         (int)aind,
                         EnvXXXCfg.ActionsV[tind][aind].dX,
                         EnvXXXCfg.ActionsV[tind][aind].dY,
                         EnvXXXCfg.ActionsV[tind][aind].endtheta,
                         sourcepose.theta * 180.0 / PI_CONST,
                         EnvXXXCfg.ActionsV[tind][aind].intermptV[EnvXXXCfg.ActionsV[tind][aind].intermptV.size() - 1].theta * 180.0 / PI_CONST,
                         EnvXXXCfg.ActionsV[tind][aind].cost,
                         mp_endx_m,
                         mp_endy_m,
                         mp_endtheta_rad);
#endif

            //add to the list of backward actions
            int targettheta = EnvXXXCfg.ActionsV[tind][aind].endtheta;
            if (targettheta < 0) targettheta = targettheta + EnvXXXCfg.NumThetaDirs;
            EnvXXXCfg.PredActionsV[targettheta].push_back(&(EnvXXXCfg.ActionsV[tind][aind]));
        }
    }

    //set number of actions
    EnvXXXCfg.actionwidth = motionprimitiveV->size();

    //now compute replanning data
    ComputeReplanningData();

    SBPL_PRINTF("done pre-computing action data based on motion primitives\n");
}

// here motionprimitivevector contains actions for all angles
void EnvironmentXXXLATTICE::PrecomputeActionswithCompleteMotionPrimitive(
    std::vector<SBPL_xxx_mprimitive>* motionprimitiveV)
{
    SBPL_PRINTF("Pre-computing action data using motion primitives for every angle...\n");
    EnvXXXCfg.ActionsV = new EnvXXXLATAction_t*[EnvXXXCfg.NumThetaDirs];
    EnvXXXCfg.PredActionsV = new std::vector<EnvXXXLATAction_t*>[EnvXXXCfg.NumThetaDirs];
    std::vector<sbpl_2Dcell_t> footprint;

    if (motionprimitiveV->size() % EnvXXXCfg.NumThetaDirs != 0) {
        throw SBPL_Exception("ERROR: motionprimitives should be uniform across actions");
    }

    EnvXXXCfg.actionwidth = ((int)motionprimitiveV->size()) / EnvXXXCfg.NumThetaDirs;

    // iterate over source angles
    int maxnumofactions = 0;
    for (int tind = 0; tind < EnvXXXCfg.NumThetaDirs; tind++) {
        SBPL_PRINTF("pre-computing for angle %d out of %d angles\n", tind, EnvXXXCfg.NumThetaDirs);

        EnvXXXCfg.ActionsV[tind] = new EnvXXXLATAction_t[EnvXXXCfg.actionwidth];

        // compute sourcepose
        sbpl_xy_theta_pt_t sourcepose;
        sourcepose.x = DISCXY2CONT(0, EnvXXXCfg.cellsize_m);
        sourcepose.y = DISCXY2CONT(0, EnvXXXCfg.cellsize_m);
        sourcepose.theta = DiscTheta2ContNew(tind);

        // iterate over motion primitives
        int numofactions = 0;
        int aind = -1;
        for (int mind = 0; mind < (int)motionprimitiveV->size(); mind++) {
            //find a motion primitive for this angle
            if (motionprimitiveV->at(mind).starttheta_c != tind) {
                continue;
            }

            aind++;
            numofactions++;

            // action index
            EnvXXXCfg.ActionsV[tind][aind].aind = aind;

            // start angle
            EnvXXXCfg.ActionsV[tind][aind].starttheta = tind;

            // compute dislocation
            EnvXXXCfg.ActionsV[tind][aind].endtheta = motionprimitiveV->at(mind).endcell.theta;
            EnvXXXCfg.ActionsV[tind][aind].dX = motionprimitiveV->at(mind).endcell.x;
            EnvXXXCfg.ActionsV[tind][aind].dY = motionprimitiveV->at(mind).endcell.y;

            // compute and store interm points as well as intersecting cells
            EnvXXXCfg.ActionsV[tind][aind].intersectingcellsV.clear();
            EnvXXXCfg.ActionsV[tind][aind].intermptV.clear();
            EnvXXXCfg.ActionsV[tind][aind].interm3DcellsV.clear();

            sbpl_xy_theta_cell_t previnterm3Dcell;
            previnterm3Dcell.x = 0;
            previnterm3Dcell.y = 0;

            // Compute all the intersected cells for this action (intermptV and interm3DcellsV)
            for (int pind = 0; pind < (int)motionprimitiveV->at(mind).intermptV.size(); pind++) {
                sbpl_xy_theta_pt_t intermpt = motionprimitiveV->at(mind).intermptV[pind];
                EnvXXXCfg.ActionsV[tind][aind].intermptV.push_back(intermpt);

                // also compute the intermediate discrete cells if not there already
                sbpl_xy_theta_pt_t pose;
                pose.x = intermpt.x + sourcepose.x;
                pose.y = intermpt.y + sourcepose.y;
                pose.theta = intermpt.theta;

                sbpl_xy_theta_cell_t intermediate2dCell;
                intermediate2dCell.x = CONTXY2DISC(pose.x, EnvXXXCfg.cellsize_m);
                intermediate2dCell.y = CONTXY2DISC(pose.y, EnvXXXCfg.cellsize_m);

                // add unique cells to the list
                if (EnvXXXCfg.ActionsV[tind][aind].interm3DcellsV.size() == 0 ||
                    intermediate2dCell.x != previnterm3Dcell.x ||
                    intermediate2dCell.y != previnterm3Dcell.y)
                {
                    EnvXXXCfg.ActionsV[tind][aind].interm3DcellsV.push_back(intermediate2dCell);
                }

                previnterm3Dcell = intermediate2dCell;
            }

            // compute linear and angular time
            double linear_distance = 0;
            for (unsigned int i = 1; i < EnvXXXCfg.ActionsV[tind][aind].intermptV.size(); i++) {
                double x0 = EnvXXXCfg.ActionsV[tind][aind].intermptV[i - 1].x;
                double y0 = EnvXXXCfg.ActionsV[tind][aind].intermptV[i - 1].y;
                double x1 = EnvXXXCfg.ActionsV[tind][aind].intermptV[i].x;
                double y1 = EnvXXXCfg.ActionsV[tind][aind].intermptV[i].y;
                double dx = x1 - x0;
                double dy = y1 - y0;
                linear_distance += sqrt(dx * dx + dy * dy);
            }
            double linear_time = linear_distance / EnvXXXCfg.nominalvel_mpersecs;
            double angular_distance;
            angular_distance = fabs(computeMinUnsignedAngleDiff(
                    DiscTheta2ContNew(EnvXXXCfg.ActionsV[tind][aind].endtheta),
                    DiscTheta2ContNew(EnvXXXCfg.ActionsV[tind][aind].starttheta)));

            double angular_time = angular_distance / ((PI_CONST / 4.0) /
                                  EnvXXXCfg.timetoturn45degsinplace_secs);
            EnvXXXCfg.ActionsV[tind][aind].time = std::max(linear_time, angular_time);
            // make the cost the max of the two times
            EnvXXXCfg.ActionsV[tind][aind].cost =
                    (int)(ceil(XXXLAT_COSTMULT_MTOMM * std::max(linear_time, angular_time)));
            // use any additional cost multiplier
            EnvXXXCfg.ActionsV[tind][aind].cost *= motionprimitiveV->at(mind).additionalactioncostmult;

            // now compute the intersecting cells for this motion (including ignoring the source footprint)
            get_2d_motion_cells(
                    EnvXXXCfg.FootprintPolygon,
                    motionprimitiveV->at(mind).intermptV,
                    &EnvXXXCfg.ActionsV[tind][aind].intersectingcellsV,
                    EnvXXXCfg.cellsize_m);

#if DEBUG
            SBPL_FPRINTF(fDeb,
                         "action tind=%2d aind=%2d: dX=%3d dY=%3d endtheta=%3d (%6.2f degs -> %6.2f degs) "
                         "cost=%4d (mprimID %3d: %3d %3d %3d) numofintermcells = %d numofintercells=%d\n",
                         tind,
                         aind,
                         EnvXXXCfg.ActionsV[tind][aind].dX,
                         EnvXXXCfg.ActionsV[tind][aind].dY,
                         EnvXXXCfg.ActionsV[tind][aind].endtheta,
                         EnvXXXCfg.ActionsV[tind][aind].intermptV[0].theta * 180 / PI_CONST,
                         EnvXXXCfg.ActionsV[tind][aind].intermptV[EnvXXXCfg.ActionsV[tind][aind].intermptV.size() - 1].theta * 180 / PI_CONST, EnvXXXCfg.ActionsV[tind][aind].cost,
                         motionprimitiveV->at(mind).motprimID, motionprimitiveV->at(mind).endcell.x,
                         motionprimitiveV->at(mind).endcell.y, motionprimitiveV->at(mind).endcell.theta,
                         (int)EnvXXXCfg.ActionsV[tind][aind].interm3DcellsV.size(),
                         (int)EnvXXXCfg.ActionsV[tind][aind].intersectingcellsV.size());
#endif

            // add to the list of backward actions
            int targettheta = EnvXXXCfg.ActionsV[tind][aind].endtheta;
            if (targettheta < 0) {
                targettheta = targettheta + EnvXXXCfg.NumThetaDirs;
            }
            EnvXXXCfg.PredActionsV[targettheta].push_back(&(EnvXXXCfg.ActionsV[tind][aind]));
        }

        if (maxnumofactions < numofactions) {
            maxnumofactions = numofactions;
        }
    }

    // at this point we don't allow nonuniform number of actions
    if (motionprimitiveV->size() != (size_t)(EnvXXXCfg.NumThetaDirs * maxnumofactions)) {
        std::stringstream ss;
        ss << "ERROR: nonuniform number of actions is not supported" <<
                " (maxnumofactions=" << maxnumofactions << " while motprims=" <<
                motionprimitiveV->size() << " thetas=" <<
                EnvXXXCfg.NumThetaDirs;
        throw SBPL_Exception(ss.str());
    }

    // now compute replanning data
    ComputeReplanningData();

    SBPL_PRINTF("done pre-computing action data based on motion primitives\n");
}

void EnvironmentXXXLATTICE::DeprecatedPrecomputeActions()
{
    SBPL_PRINTF("Use of DeprecatedPrecomputeActions() is deprecated and probably doesn't work!\n");

    // construct list of actions
    SBPL_PRINTF("Pre-computing action data using the motion primitives for a 3D kinematic planning...\n");
    EnvXXXCfg.ActionsV = new EnvXXXLATAction_t*[EnvXXXCfg.NumThetaDirs];
    EnvXXXCfg.PredActionsV = new std::vector<EnvXXXLATAction_t*> [EnvXXXCfg.NumThetaDirs];
    std::vector<sbpl_2Dcell_t> footprint;
    // iterate over source angles
    for (int tind = 0; tind < EnvXXXCfg.NumThetaDirs; tind++) {
        SBPL_PRINTF("processing angle %d\n", tind);
        EnvXXXCfg.ActionsV[tind] = new EnvXXXLATAction_t[EnvXXXCfg.actionwidth];

        // compute sourcepose
        sbpl_xy_theta_pt_t sourcepose;
        sourcepose.x = DISCXY2CONT(0, EnvXXXCfg.cellsize_m);
        sourcepose.y = DISCXY2CONT(0, EnvXXXCfg.cellsize_m);
        sourcepose.theta = DiscTheta2ContNew(tind);

        // the construction assumes that the robot first turns and then goes
        // along this new theta
        int aind = 0;
        for (; aind < 3; aind++) {
            EnvXXXCfg.ActionsV[tind][aind].aind = aind;
            EnvXXXCfg.ActionsV[tind][aind].starttheta = tind;
            // -1,0,1
            EnvXXXCfg.ActionsV[tind][aind].endtheta = (tind + aind - 1) % EnvXXXCfg.NumThetaDirs;
            double angle = DiscTheta2ContNew(EnvXXXCfg.ActionsV[tind][aind].endtheta);
            EnvXXXCfg.ActionsV[tind][aind].dX = (int)(cos(angle) + 0.5 * (cos(angle) > 0 ? 1 : -1));
            EnvXXXCfg.ActionsV[tind][aind].dY = (int)(sin(angle) + 0.5 * (sin(angle) > 0 ? 1 : -1));
            EnvXXXCfg.ActionsV[tind][aind].cost = (int)(ceil(XXXLAT_COSTMULT_MTOMM *
                    EnvXXXCfg.cellsize_m / EnvXXXCfg.nominalvel_mpersecs *
                    sqrt((double)(EnvXXXCfg.ActionsV[tind][aind].dX *
                            EnvXXXCfg.ActionsV[tind][aind].dX + EnvXXXCfg.ActionsV[tind][aind].dY *
                            EnvXXXCfg.ActionsV[tind][aind].dY))));

            // compute intersecting cells
            sbpl_xy_theta_pt_t pose;
            pose.x = DISCXY2CONT(EnvXXXCfg.ActionsV[tind][aind].dX, EnvXXXCfg.cellsize_m);
            pose.y = DISCXY2CONT(EnvXXXCfg.ActionsV[tind][aind].dY, EnvXXXCfg.cellsize_m);
            pose.theta = angle;
            EnvXXXCfg.ActionsV[tind][aind].intermptV.clear();
            EnvXXXCfg.ActionsV[tind][aind].intersectingcellsV.clear();
            get_2d_footprint_cells(
                    EnvXXXCfg.FootprintPolygon,
                    &EnvXXXCfg.ActionsV[tind][aind].intersectingcellsV,
                    pose,
                    EnvXXXCfg.cellsize_m);
            RemoveSourceFootprint(sourcepose, &EnvXXXCfg.ActionsV[tind][aind].intersectingcellsV);

#if DEBUG
            SBPL_PRINTF("action tind=%d aind=%d: endtheta=%d (%f) dX=%d dY=%d cost=%d\n",
                        tind, aind, EnvXXXCfg.ActionsV[tind][aind].endtheta, angle,
                        EnvXXXCfg.ActionsV[tind][aind].dX, EnvXXXCfg.ActionsV[tind][aind].dY,
                        EnvXXXCfg.ActionsV[tind][aind].cost);
#endif

            // add to the list of backward actions
            int targettheta = EnvXXXCfg.ActionsV[tind][aind].endtheta;
            if (targettheta < 0) {
                targettheta = targettheta + EnvXXXCfg.NumThetaDirs;
            }
            EnvXXXCfg.PredActionsV[targettheta].push_back(&(EnvXXXCfg.ActionsV[tind][aind]));
        }

        // decrease and increase angle without movement
        aind = 3;
        EnvXXXCfg.ActionsV[tind][aind].aind = aind;
        EnvXXXCfg.ActionsV[tind][aind].starttheta = tind;
        EnvXXXCfg.ActionsV[tind][aind].endtheta = tind - 1;
        if (EnvXXXCfg.ActionsV[tind][aind].endtheta < 0) {
            EnvXXXCfg.ActionsV[tind][aind].endtheta += EnvXXXCfg.NumThetaDirs;
        }
        EnvXXXCfg.ActionsV[tind][aind].dX = 0;
        EnvXXXCfg.ActionsV[tind][aind].dY = 0;
        EnvXXXCfg.ActionsV[tind][aind].cost =
                (int)(XXXLAT_COSTMULT_MTOMM * EnvXXXCfg.timetoturn45degsinplace_secs);

        // compute intersecting cells
        sbpl_xy_theta_pt_t pose;
        pose.x = DISCXY2CONT(EnvXXXCfg.ActionsV[tind][aind].dX, EnvXXXCfg.cellsize_m);
        pose.y = DISCXY2CONT(EnvXXXCfg.ActionsV[tind][aind].dY, EnvXXXCfg.cellsize_m);
        pose.theta = DiscTheta2ContNew(EnvXXXCfg.ActionsV[tind][aind].endtheta);
        EnvXXXCfg.ActionsV[tind][aind].intermptV.clear();
        EnvXXXCfg.ActionsV[tind][aind].intersectingcellsV.clear();
        get_2d_footprint_cells(
                EnvXXXCfg.FootprintPolygon,
                &EnvXXXCfg.ActionsV[tind][aind].intersectingcellsV,
                pose,
                EnvXXXCfg.cellsize_m);
        RemoveSourceFootprint(sourcepose, &EnvXXXCfg.ActionsV[tind][aind].intersectingcellsV);

#if DEBUG
        SBPL_PRINTF("action tind=%d aind=%d: endtheta=%d (%f) dX=%d dY=%d cost=%d\n",
                    tind, aind, EnvXXXCfg.ActionsV[tind][aind].endtheta,
                    DiscTheta2ContNew(EnvXXXCfg.ActionsV[tind][aind].endtheta),
                    EnvXXXCfg.ActionsV[tind][aind].dX, EnvXXXCfg.ActionsV[tind][aind].dY,
                    EnvXXXCfg.ActionsV[tind][aind].cost);
#endif

        // add to the list of backward actions
        int targettheta = EnvXXXCfg.ActionsV[tind][aind].endtheta;
        if (targettheta < 0) {
            targettheta = targettheta + EnvXXXCfg.NumThetaDirs;
        }
        EnvXXXCfg.PredActionsV[targettheta].push_back(&(EnvXXXCfg.ActionsV[tind][aind]));

        aind = 4;
        EnvXXXCfg.ActionsV[tind][aind].aind = aind;
        EnvXXXCfg.ActionsV[tind][aind].starttheta = tind;
        EnvXXXCfg.ActionsV[tind][aind].endtheta = (tind + 1) % EnvXXXCfg.NumThetaDirs;
        EnvXXXCfg.ActionsV[tind][aind].dX = 0;
        EnvXXXCfg.ActionsV[tind][aind].dY = 0;
        EnvXXXCfg.ActionsV[tind][aind].cost =
                (int)(XXXLAT_COSTMULT_MTOMM * EnvXXXCfg.timetoturn45degsinplace_secs);

        // compute intersecting cells
        pose.x = DISCXY2CONT(EnvXXXCfg.ActionsV[tind][aind].dX, EnvXXXCfg.cellsize_m);
        pose.y = DISCXY2CONT(EnvXXXCfg.ActionsV[tind][aind].dY, EnvXXXCfg.cellsize_m);
        pose.theta = DiscTheta2ContNew(EnvXXXCfg.ActionsV[tind][aind].endtheta);
        EnvXXXCfg.ActionsV[tind][aind].intermptV.clear();
        EnvXXXCfg.ActionsV[tind][aind].intersectingcellsV.clear();
        get_2d_footprint_cells(
                EnvXXXCfg.FootprintPolygon,
                &EnvXXXCfg.ActionsV[tind][aind].intersectingcellsV,
                pose,
                EnvXXXCfg.cellsize_m);
        RemoveSourceFootprint(sourcepose, &EnvXXXCfg.ActionsV[tind][aind].intersectingcellsV);

#if DEBUG
        SBPL_PRINTF("action tind=%d aind=%d: endtheta=%d (%f) dX=%d dY=%d cost=%d\n",
                    tind, aind, EnvXXXCfg.ActionsV[tind][aind].endtheta, DiscTheta2ContNew(EnvXXXCfg.ActionsV[tind][aind].endtheta),
                    EnvXXXCfg.ActionsV[tind][aind].dX, EnvXXXCfg.ActionsV[tind][aind].dY,
                    EnvXXXCfg.ActionsV[tind][aind].cost);
#endif

        // add to the list of backward actions
        targettheta = EnvXXXCfg.ActionsV[tind][aind].endtheta;
        if (targettheta < 0) targettheta = targettheta + EnvXXXCfg.NumThetaDirs;
        EnvXXXCfg.PredActionsV[targettheta].push_back(&(EnvXXXCfg.ActionsV[tind][aind]));
    }

    // now compute replanning data
    ComputeReplanningData();

    SBPL_PRINTF("done pre-computing action data\n");
}

void EnvironmentXXXLATTICE::InitializeEnvConfig(std::vector<SBPL_xxx_mprimitive>* motionprimitiveV)
{
    // additional to configuration file initialization of EnvXXXCfg if
    // necessary

    // dXY dirs
    EnvXXXCfg.dXY[0][0] = -1;
    EnvXXXCfg.dXY[0][1] = -1;
    EnvXXXCfg.dXY[1][0] = -1;
    EnvXXXCfg.dXY[1][1] = 0;
    EnvXXXCfg.dXY[2][0] = -1;
    EnvXXXCfg.dXY[2][1] = 1;
    EnvXXXCfg.dXY[3][0] = 0;
    EnvXXXCfg.dXY[3][1] = -1;
    EnvXXXCfg.dXY[4][0] = 0;
    EnvXXXCfg.dXY[4][1] = 1;
    EnvXXXCfg.dXY[5][0] = 1;
    EnvXXXCfg.dXY[5][1] = -1;
    EnvXXXCfg.dXY[6][0] = 1;
    EnvXXXCfg.dXY[6][1] = 0;
    EnvXXXCfg.dXY[7][0] = 1;
    EnvXXXCfg.dXY[7][1] = 1;

    sbpl_xy_theta_pt_t temppose;
    temppose.x = 0.0;
    temppose.y = 0.0;
    temppose.theta = 0.0;
    std::vector<sbpl_2Dcell_t> footprint;
    get_2d_footprint_cells(
            EnvXXXCfg.FootprintPolygon,
            &footprint,
            temppose,
            EnvXXXCfg.cellsize_m);
    SBPL_PRINTF("number of cells in footprint of the robot = %d\n", (unsigned int)footprint.size());

    // for (std::vector<sbpl_2Dcell_t>::iterator it = footprint.begin(); it != footprint.end(); ++it) {
    //     SBPL_PRINTF("Footprint cell at (%d, %d)\n", it->x, it->y);
    // }

#if DEBUG
    SBPL_FPRINTF(fDeb, "footprint cells (size=%d):\n", (int)footprint.size());
    for(int i = 0; i < (int) footprint.size(); i++)
    {
        SBPL_FPRINTF(fDeb, "%d %d (cont: %.3f %.3f)\n", footprint.at(i).x, footprint.at(i).y,
                     DISCXY2CONT(footprint.at(i).x, EnvXXXCfg.cellsize_m),
                     DISCXY2CONT(footprint.at(i).y, EnvXXXCfg.cellsize_m));
    }
#endif

    if (motionprimitiveV == NULL) {
        DeprecatedPrecomputeActions();
    }
    else {
        PrecomputeActionswithCompleteMotionPrimitive(motionprimitiveV);
    }
}

bool EnvironmentXXXLATTICE::IsValidCell(int X, int Y)
{
    return (X >= 0 && X < EnvXXXCfg.EnvWidth_c &&
            Y >= 0 && Y < EnvXXXCfg.EnvHeight_c &&
            EnvXXXCfg.Grid2D[X][Y] < EnvXXXCfg.obsthresh);
}

bool EnvironmentXXXLATTICE::IsWithinMapCell(int X, int Y)
{
    return (X >= 0 && X < EnvXXXCfg.EnvWidth_c &&
            Y >= 0 && Y < EnvXXXCfg.EnvHeight_c);
}

bool EnvironmentXXXLATTICE::IsValidConfiguration(int X, int Y, int Theta)
{
    std::vector<sbpl_2Dcell_t> footprint;
    sbpl_xy_theta_pt_t pose;

    // compute continuous pose
    pose.x = DISCXY2CONT(X, EnvXXXCfg.cellsize_m);
    pose.y = DISCXY2CONT(Y, EnvXXXCfg.cellsize_m);
    pose.theta = DiscTheta2ContNew(Theta);

    // compute footprint cells
    get_2d_footprint_cells(
            EnvXXXCfg.FootprintPolygon,
            &footprint,
            pose,
            EnvXXXCfg.cellsize_m);

    // iterate over all footprint cells
    for (int find = 0; find < (int)footprint.size(); find++) {
        int x = footprint.at(find).x;
        int y = footprint.at(find).y;

        if (x < 0 || x >= EnvXXXCfg.EnvWidth_c ||
            y < 0 || y >= EnvXXXCfg.EnvHeight_c ||
            EnvXXXCfg.Grid2D[x][y] >= EnvXXXCfg.obsthresh)
        {
            return false;
        }
    }

    return true;
}

int EnvironmentXXXLATTICE::GetActionCost(
    int SourceX, int SourceY, int SourceTheta,
    EnvXXXLATAction_t* action)
{
    sbpl_2Dcell_t cell;
    sbpl_xy_theta_cell_t interm3Dcell;
    int i;

    // TODO - go over bounding box (minpt and maxpt) to test validity and skip
    // testing boundaries below, also order intersect cells so that the four
    // farthest pts go first

    if (!IsValidCell(SourceX, SourceY)) {
        return INFINITECOST;
    }
    if (!IsValidCell(SourceX + action->dX, SourceY + action->dY)) {
        return INFINITECOST;
    }

    if (EnvXXXCfg.Grid2D[SourceX + action->dX][SourceY + action->dY] >=
        EnvXXXCfg.cost_inscribed_thresh)
    {
        return INFINITECOST;
    }

    // need to iterate over discretized center cells and compute cost based on them
    unsigned char maxcellcost = 0;
    for (i = 0; i < (int)action->interm3DcellsV.size(); i++) {
        interm3Dcell = action->interm3DcellsV.at(i);
        interm3Dcell.x = interm3Dcell.x + SourceX;
        interm3Dcell.y = interm3Dcell.y + SourceY;

        if (interm3Dcell.x < 0 || interm3Dcell.x >= EnvXXXCfg.EnvWidth_c ||
            interm3Dcell.y < 0 || interm3Dcell.y >= EnvXXXCfg.EnvHeight_c)
        {
            return INFINITECOST;
        }

        maxcellcost = __max(maxcellcost, EnvXXXCfg.Grid2D[interm3Dcell.x][interm3Dcell.y]);

        // check that the robot is NOT in the cell at which there is no valid orientation
        if (maxcellcost >= EnvXXXCfg.cost_inscribed_thresh) {
            return INFINITECOST;
        }
    }

    // check collisions that for the particular footprint orientation along the action
    if (EnvXXXCfg.FootprintPolygon.size() > 1 &&
        (int)maxcellcost >= EnvXXXCfg.cost_possibly_circumscribed_thresh)
    {
        checks++;

        for (i = 0; i < (int)action->intersectingcellsV.size(); i++) {
            // get the cell in the map
            cell = action->intersectingcellsV.at(i);
            cell.x = cell.x + SourceX;
            cell.y = cell.y + SourceY;

            // check validity
            if (!IsValidCell(cell.x, cell.y)) {
                return INFINITECOST;
            }

// cost computation changed: cost = max(cost of centers of the robot along
// action) intersecting cells are only used for collision checking
//            if (EnvXXXCfg.Grid2D[cell.x][cell.y] > currentmaxcost) {
//              currentmaxcost = EnvXXXCfg.Grid2D[cell.x][cell.y];
//            }
        }
    }

    // to ensure consistency of h2D:
    maxcellcost = __max(maxcellcost, EnvXXXCfg.Grid2D[SourceX][SourceY]);
    int currentmaxcost = (int)__max(
            maxcellcost,
            EnvXXXCfg.Grid2D[SourceX + action->dX][SourceY + action->dY]);

    // use cell cost as multiplicative factor
    return action->cost * (currentmaxcost + 1);
}

double EnvironmentXXXLATTICE::EuclideanDistance_m(int X1, int Y1, int X2, int Y2)
{
    int sqdist = ((X1 - X2) * (X1 - X2) + (Y1 - Y2) * (Y1 - Y2));
    return EnvXXXCfg.cellsize_m * sqrt((double)sqdist);
}

// calculates a set of cells that correspond to the specified footprint adds
// points to it (does not clear it beforehand)
void EnvironmentXXXLATTICE::CalculateFootprintForPose(
    sbpl_xy_theta_pt_t pose, std::vector<sbpl_2Dcell_t>* footprint,
    const std::vector<sbpl_2Dpt_t>& FootprintPolygon)
{
    int pind;

    // handle special case where footprint is just a point
    if (FootprintPolygon.size() <= 1) {
        sbpl_2Dcell_t cell;
        cell.x = CONTXY2DISC(pose.x, EnvXXXCfg.cellsize_m);
        cell.y = CONTXY2DISC(pose.y, EnvXXXCfg.cellsize_m);

        for (pind = 0; pind < (int)footprint->size(); pind++) {
            if (cell.x == footprint->at(pind).x && cell.y == footprint->at(pind).y) break;
        }
        if (pind == (int)footprint->size()) {
            footprint->push_back(cell);
        }
        return;
    }

    std::vector<sbpl_2Dpt_t> bounding_polygon;
    unsigned int find;
    double max_x = -INFINITECOST, min_x = INFINITECOST, max_y = -INFINITECOST, min_y = INFINITECOST;
    sbpl_2Dpt_t pt(0, 0);
    for (find = 0; find < FootprintPolygon.size(); find++) {
        // rotate and translate the corner of the robot
        pt = FootprintPolygon[find];

        // rotate and translate the point
        sbpl_2Dpt_t corner;
        corner.x = cos(pose.theta) * pt.x - sin(pose.theta) * pt.y + pose.x;
        corner.y = sin(pose.theta) * pt.x + cos(pose.theta) * pt.y + pose.y;
        bounding_polygon.push_back(corner);
        if (corner.x < min_x || find == 0) {
            min_x = corner.x;
        }
        if (corner.x > max_x || find == 0) {
            max_x = corner.x;
        }
        if (corner.y < min_y || find == 0) {
            min_y = corner.y;
        }
        if (corner.y > max_y || find == 0) {
            max_y = corner.y;
        }
    }

    // initialize previous values to something that will fail the if condition
    // during the first iteration in the for loop
    int prev_discrete_x = CONTXY2DISC(pt.x, EnvXXXCfg.cellsize_m) + 1;
    int prev_discrete_y = CONTXY2DISC(pt.y, EnvXXXCfg.cellsize_m) + 1;
    int prev_inside = 0;
    int discrete_x;
    int discrete_y;

    for (double x = min_x; x <= max_x; x += EnvXXXCfg.cellsize_m / 3) {
        for (double y = min_y; y <= max_y; y += EnvXXXCfg.cellsize_m / 3) {
            pt.x = x;
            pt.y = y;
            discrete_x = CONTXY2DISC(pt.x, EnvXXXCfg.cellsize_m);
            discrete_y = CONTXY2DISC(pt.y, EnvXXXCfg.cellsize_m);

            // see if we just tested this point
            if (discrete_x != prev_discrete_x || discrete_y != prev_discrete_y || prev_inside == 0) {

                if (IsInsideFootprint(pt, &bounding_polygon)) {
                    // convert to a grid point

                    sbpl_2Dcell_t cell;
                    cell.x = discrete_x;
                    cell.y = discrete_y;

                    // insert point if not there already
                    for (pind = 0; pind < (int)footprint->size(); pind++) {
                        if (cell.x == footprint->at(pind).x && cell.y == footprint->at(pind).y) break;
                    }
                    if (pind == (int)footprint->size()) {
                        footprint->push_back(cell);
                    }

                    prev_inside = 1;

                }
                else {
                    prev_inside = 0;
                }

            }
            else {

            }

            prev_discrete_x = discrete_x;
            prev_discrete_y = discrete_y;
        } // over x_min...x_max
    }
}

// calculates a set of cells that correspond to the footprint of the base adds
// points to it (does not clear it beforehand)
void EnvironmentXXXLATTICE::CalculateFootprintForPose(
    sbpl_xy_theta_pt_t pose,
    std::vector<sbpl_2Dcell_t>* footprint)
{
    CalculateFootprintForPose(pose, footprint, EnvXXXCfg.FootprintPolygon);
}

// removes a set of cells that correspond to the specified footprint at the
// sourcepose adds points to it (does not clear it beforehand)
void EnvironmentXXXLATTICE::RemoveSourceFootprint(
    sbpl_xy_theta_pt_t sourcepose,
    std::vector<sbpl_2Dcell_t>* footprint,
    const std::vector<sbpl_2Dpt_t>& FootprintPolygon)
{
    std::vector<sbpl_2Dcell_t> sourcefootprint;

    // compute source footprint
    get_2d_footprint_cells(FootprintPolygon, &sourcefootprint, sourcepose, EnvXXXCfg.cellsize_m);

    // now remove the source cells from the footprint
    for (int sind = 0; sind < (int)sourcefootprint.size(); sind++) {
        for (int find = 0; find < (int)footprint->size(); find++) {
            if (sourcefootprint.at(sind).x == footprint->at(find).x && sourcefootprint.at(sind).y
                == footprint->at(find).y) {
                footprint->erase(footprint->begin() + find);
                break;
            }
        } // over footprint
    } // over source
}

// removes a set of cells that correspond to the footprint of the base at the
// sourcepose adds points to it (does not clear it beforehand)
void EnvironmentXXXLATTICE::RemoveSourceFootprint(
    sbpl_xy_theta_pt_t sourcepose,
    std::vector<sbpl_2Dcell_t>* footprint)
{
    RemoveSourceFootprint(sourcepose, footprint, EnvXXXCfg.FootprintPolygon);
}

void EnvironmentXXXLATTICE::EnsureHeuristicsUpdated(bool bGoalHeuristics)
{
    if (bNeedtoRecomputeStartHeuristics && !bGoalHeuristics) {
        grid2Dsearchfromstart->search(
                EnvXXXCfg.Grid2D,
                EnvXXXCfg.cost_inscribed_thresh,
                EnvXXXCfg.StartX_c, EnvXXXCfg.StartY_c,
                EnvXXXCfg.EndX_c, EnvXXXCfg.EndY_c,
                SBPL_2DGRIDSEARCH_TERM_CONDITION_TWOTIMESOPTPATH);
        bNeedtoRecomputeStartHeuristics = false;
        SBPL_PRINTF("2dsolcost_infullunits=%d\n", (int)(grid2Dsearchfromstart->getlowerboundoncostfromstart_inmm(EnvXXXCfg.EndX_c, EnvXXXCfg.EndY_c) / EnvXXXCfg.nominalvel_mpersecs));
    }

    if (bNeedtoRecomputeGoalHeuristics && bGoalHeuristics) {
        grid2Dsearchfromgoal->search(
                EnvXXXCfg.Grid2D,
                EnvXXXCfg.cost_inscribed_thresh,
                EnvXXXCfg.EndX_c, EnvXXXCfg.EndY_c,
                EnvXXXCfg.StartX_c, EnvXXXCfg.StartY_c,
                SBPL_2DGRIDSEARCH_TERM_CONDITION_TWOTIMESOPTPATH);
        bNeedtoRecomputeGoalHeuristics = false;
        SBPL_PRINTF("2dsolcost_infullunits=%d\n", (int)(grid2Dsearchfromgoal->getlowerboundoncostfromstart_inmm(EnvXXXCfg.StartX_c, EnvXXXCfg.StartY_c) / EnvXXXCfg.nominalvel_mpersecs));
    }
}

void EnvironmentXXXLATTICE::ComputeHeuristicValues()
{
    // whatever necessary pre-computation of heuristic values is done here
    SBPL_PRINTF("Precomputing heuristics...\n");

    // allocated 2D grid searches
    grid2Dsearchfromstart = new SBPL2DGridSearch(
            EnvXXXCfg.EnvWidth_c, EnvXXXCfg.EnvHeight_c,
            (float)EnvXXXCfg.cellsize_m, blocksize, bucketsize);
    grid2Dsearchfromgoal = new SBPL2DGridSearch(
            EnvXXXCfg.EnvWidth_c, EnvXXXCfg.EnvHeight_c,
            (float)EnvXXXCfg.cellsize_m, blocksize, bucketsize);

    // set OPEN type to sliding buckets
    grid2Dsearchfromstart->setOPENdatastructure(SBPL_2DGRIDSEARCH_OPENTYPE_SLIDINGBUCKETS);
    grid2Dsearchfromgoal->setOPENdatastructure(SBPL_2DGRIDSEARCH_OPENTYPE_SLIDINGBUCKETS);

    SBPL_PRINTF("done\n");
}

bool EnvironmentXXXLATTICE::CheckQuant(FILE* fOut)
{
    for (double theta = -10; theta < 10;
        theta += 2.0 * PI_CONST / EnvXXXCfg.NumThetaDirs * 0.01)
    {
        int nTheta, nnewTheta;
        double newTheta;
        nTheta = ContTheta2DiscNew(theta);
        newTheta = DiscTheta2ContNew(nTheta);
        nnewTheta = ContTheta2DiscNew(newTheta);

        SBPL_FPRINTF(fOut, "theta=%f(%f)->%d->%f->%d\n", theta, theta * 180 / PI_CONST, nTheta, newTheta, nnewTheta);

        if (nTheta != nnewTheta) {
            SBPL_ERROR("ERROR: invalid quantization\n");
            return false;
        }
    }

    return true;
}

bool EnvironmentXXXLATTICE::InitializeEnv(
    const char* sEnvFile,
    const std::vector<sbpl_2Dpt_t>& perimeterptsV,
    const std::vector<sbpl_2Dpt_t>& trailer_perimeterptsV,
    const char* sMotPrimFile)
{
    EnvXXXCfg.FootprintPolygon = perimeterptsV;
    EnvXXXCfg.TrailerPolygon = trailer_perimeterptsV;

    SBPL_INFO("InitializeEnv start: sEnvFile=%s sMotPrimFile=%s\n", sEnvFile, sMotPrimFile);
    fflush(stdout);

    FILE* fCfg = fopen(sEnvFile, "r");
    if (fCfg == NULL) {
        std::stringstream ss;
        ss << "ERROR: unable to open " << sEnvFile;
        throw SBPL_Exception(ss.str());
    }

    ReadConfiguration(fCfg);
    fclose(fCfg);

    if (sMotPrimFile != NULL) {
        FILE* fMotPrim = fopen(sMotPrimFile, "r");
        if (fMotPrim == NULL) {
            std::stringstream ss;
            ss << "ERROR: unable to open " << sMotPrimFile;
            throw SBPL_Exception(ss.str());
        }
        if (ReadMotionPrimitives(fMotPrim) == false) {
            throw SBPL_Exception("ERROR: failed to read in motion primitive file");
        }

        EnvXXXCfg.StartTheta = ContTheta2DiscNew(EnvXXXCfg.StartTheta_rad);
        if (EnvXXXCfg.StartTheta < 0 ||
            EnvXXXCfg.StartTheta >= EnvXXXCfg.NumThetaDirs)
        {
            throw new SBPL_Exception("ERROR: illegal start coordinates for theta");
        }
        EnvXXXCfg.EndTheta = ContTheta2DiscNew(EnvXXXCfg.EndTheta_rad);
        if (EnvXXXCfg.EndTheta < 0 ||
            EnvXXXCfg.EndTheta >= EnvXXXCfg.NumThetaDirs)
        {
            throw new SBPL_Exception("ERROR: illegal goal coordinates for theta");
        }

        InitGeneral(&EnvXXXCfg.mprimV);
        fclose(fMotPrim);
    }
    else {
        InitGeneral( NULL);
    }

    SBPL_PRINTF("size of env: %d by %d\n", EnvXXXCfg.EnvWidth_c, EnvXXXCfg.EnvHeight_c);
    // stateParents.clear();
    // trailerTransitions.clear();

    return true;
}

bool EnvironmentXXXLATTICE::InitializeEnv(const char* sEnvFile)
{
    FILE* fCfg = fopen(sEnvFile, "r");
    if (fCfg == NULL) {
        SBPL_ERROR("ERROR: unable to open %s\n", sEnvFile);
        throw SBPL_Exception();
    }
    ReadConfiguration(fCfg);
    fclose(fCfg);

    InitGeneral( NULL);

    return true;
}

bool EnvironmentXXXLATTICE::InitializeEnv(
    int width,
    int height,
    const std::vector<sbpl_2Dpt_t>& perimeterptsV,
    const std::vector<sbpl_2Dpt_t>& trailer_perimeterptsV,
    double cellsize_m,
    double nominalvel_mpersecs,
    double timetoturn45degsinplace_secs,
    unsigned char obsthresh,
    const char* sMotPrimFile,
    EnvXXXLAT_InitParms params)
{
    EnvXXXCfg.NumThetaDirs = params.numThetas;

    return InitializeEnv(
            width, height,
            params.mapdata,
            params.startx, params.starty, params.starttheta, params.starttheta1, params.starttheta2,
            params.goalx, params.goaly, params.goaltheta,
            params.goaltol_x, params.goaltol_y, params.goaltol_theta,
            perimeterptsV,
            trailer_perimeterptsV,
            cellsize_m,
            nominalvel_mpersecs, timetoturn45degsinplace_secs,
            obsthresh,
            sMotPrimFile);
}

bool EnvironmentXXXLATTICE::InitializeEnv(
    int width,
    int height,
    const unsigned char* mapdata,
    double startx, double starty, double starttheta, double starttheta1, double starttheta2,
    double goalx, double goaly, double goaltheta,
    double goaltol_x, double goaltol_y, double goaltol_theta,
    const std::vector<sbpl_2Dpt_t> & perimeterptsV,
    const std::vector<sbpl_2Dpt_t> & trailer_perimeterptsV,
    double cellsize_m,
    double nominalvel_mpersecs,
    double timetoturn45degsinplace_secs,
    unsigned char obsthresh,
    const char* sMotPrimFile)
{
    SBPL_PRINTF("env: initialize with width=%d height=%d start=%.3f %.3f %.3f "
                "goalx=%.3f %.3f %.3f cellsize=%.3f nomvel=%.3f timetoturn=%.3f, obsthresh=%d\n",
                width, height, startx, starty, starttheta, goalx, goaly, goaltheta, cellsize_m, nominalvel_mpersecs,
                timetoturn45degsinplace_secs, obsthresh);

    SBPL_PRINTF("NOTE: goaltol parameters currently unused\n");

    SBPL_PRINTF("perimeter has size=%d\n", (unsigned int)perimeterptsV.size());

    for (int i = 0; i < (int)perimeterptsV.size(); i++) {
        SBPL_PRINTF("perimeter(%d) = %.4f %.4f\n", i, perimeterptsV.at(i).x, perimeterptsV.at(i).y);
    }

    EnvXXXCfg.obsthresh = obsthresh;
    EnvXXXCfg.cellsize_m = cellsize_m;
    EnvXXXCfg.StartTheta_rad = starttheta;
    EnvXXXCfg.StartTheta1 = starttheta1;
    EnvXXXCfg.StartTheta2 = starttheta2;
    EnvXXXCfg.EndTheta_rad = goaltheta;

    // TODO - need to set the tolerance as well

    if (sMotPrimFile != NULL) {
        FILE* fMotPrim = fopen(sMotPrimFile, "r");
        if (fMotPrim == NULL) {
            std::stringstream ss;
            ss << "ERROR: unable to open " << sMotPrimFile;
            throw SBPL_Exception(ss.str());
        }

        if (ReadMotionPrimitives(fMotPrim) == false) {
            throw SBPL_Exception("ERROR: failed to read in motion primitive file");
        }
        fclose(fMotPrim);
    }

    EnvXXXCfg.StartTheta = ContTheta2DiscNew(EnvXXXCfg.StartTheta_rad);
    if (EnvXXXCfg.StartTheta < 0 ||
        EnvXXXCfg.StartTheta >= EnvXXXCfg.NumThetaDirs)
    {
        throw new SBPL_Exception("ERROR: illegal start coordinates for theta");
    }
    EnvXXXCfg.EndTheta = ContTheta2DiscNew(EnvXXXCfg.EndTheta_rad);
    if (EnvXXXCfg.EndTheta < 0 ||
        EnvXXXCfg.EndTheta >= EnvXXXCfg.NumThetaDirs)
    {
        throw new SBPL_Exception("ERROR: illegal goal coordiantes for theta");
    }

    SetConfiguration(
            width, height, mapdata,
            CONTXY2DISC(startx, cellsize_m), CONTXY2DISC(starty, cellsize_m), EnvXXXCfg.StartTheta, starttheta1, starttheta2,
            CONTXY2DISC(goalx, cellsize_m), CONTXY2DISC(goaly, cellsize_m), EnvXXXCfg.EndTheta,
            cellsize_m,
            nominalvel_mpersecs, timetoturn45degsinplace_secs,
            perimeterptsV,
            trailer_perimeterptsV);

    if (EnvXXXCfg.mprimV.size() != 0) {
        InitGeneral(&EnvXXXCfg.mprimV);
    }
    else {
        InitGeneral( NULL);
    }

    return true;
}

bool EnvironmentXXXLATTICE::InitGeneral(
    std::vector<SBPL_xxx_mprimitive>* motionprimitiveV)
{
    // Initialize other parameters of the environment
    InitializeEnvConfig(motionprimitiveV);

    // initialize Environment
    InitializeEnvironment();

    // pre-compute heuristics
    ComputeHeuristicValues();

    return true;
}

bool EnvironmentXXXLATTICE::InitializeMDPCfg(MDPConfig* MDPCfg)
{
    // initialize MDPCfg with the start and goal ids
    MDPCfg->goalstateid = EnvXXXLAT.goalstateid;
    MDPCfg->startstateid = EnvXXXLAT.startstateid;

    return true;
}

void EnvironmentXXXLATTICE::PrintHeuristicValues()
{
    const char* heur = "heur.txt";
    FILE* fHeur = SBPL_FOPEN(heur, "w");
    if (fHeur == NULL) {
        throw SBPL_Exception("ERROR: could not open debug file to write heuristic");
    }
    SBPL2DGridSearch* grid2Dsearch = NULL;

    for (int i = 0; i < 2; i++) {
        if (i == 0 && grid2Dsearchfromstart != NULL) {
            grid2Dsearch = grid2Dsearchfromstart;
            SBPL_FPRINTF(fHeur, "start heuristics:\n");
        }
        else if (i == 1 && grid2Dsearchfromgoal != NULL) {
            grid2Dsearch = grid2Dsearchfromgoal;
            SBPL_FPRINTF(fHeur, "goal heuristics:\n");
        }
        else {
            continue;
        }

        for (int y = 0; y < EnvXXXCfg.EnvHeight_c; y++) {
            for (int x = 0; x < EnvXXXCfg.EnvWidth_c; x++) {
                if (grid2Dsearch->getlowerboundoncostfromstart_inmm(x, y) < INFINITECOST) {
                    SBPL_FPRINTF(fHeur, "%5d ", grid2Dsearch->getlowerboundoncostfromstart_inmm(x, y));
                }
                else {
                    SBPL_FPRINTF(fHeur, "XXXXX ");
                }
            }
            SBPL_FPRINTF(fHeur, "\n");
        }
    }
    SBPL_FCLOSE(fHeur);
}

void EnvironmentXXXLATTICE::SetAllPreds(CMDPSTATE* state)
{
    // implement this if the planner needs access to predecessors
    throw SBPL_Exception("ERROR in EnvXXXLAT... function: SetAllPreds is undefined");
}

void EnvironmentXXXLATTICE::GetSuccs(
    int SourceStateID,
    std::vector<int>* SuccIDV,
    std::vector<int>* CostV)
{
    GetSuccs(SourceStateID, SuccIDV, CostV, NULL);
}
void EnvironmentXXXLATTICE::GetLazySuccs(
    int SourceStateID,
    std::vector<int>* SuccIDV,
    std::vector<int>* CostV,
    std::vector<bool>* isTrueCost)
{
    GetLazySuccs(SourceStateID, SuccIDV, CostV, isTrueCost, NULL);
}

void EnvironmentXXXLATTICE::GetSuccsWithUniqueIds(
    int SourceStateID,
    std::vector<int>* SuccIDV,
    std::vector<int>* CostV)
{
    GetSuccsWithUniqueIds(SourceStateID, SuccIDV, CostV, NULL);
}

void EnvironmentXXXLATTICE::GetLazySuccsWithUniqueIds(
    int SourceStateID,
    std::vector<int>* SuccIDV,
    std::vector<int>* CostV,
    std::vector<bool>* isTrueCost)
{
    GetLazySuccsWithUniqueIds(SourceStateID, SuccIDV, CostV, isTrueCost, NULL);
}

const EnvXXXLATConfig_t* EnvironmentXXXLATTICE::GetEnvNavConfig()
{
    return &EnvXXXCfg;
}

bool EnvironmentXXXLATTICE::UpdateCost(
    int x,
    int y,
    unsigned char newcost)
{
    EnvXXXCfg.Grid2D[x][y] = newcost;

    bNeedtoRecomputeStartHeuristics = true;
    bNeedtoRecomputeGoalHeuristics = true;

    return true;
}

bool EnvironmentXXXLATTICE::SetMap(const unsigned char* mapdata)
{
    int xind = -1, yind = -1;

    for (xind = 0; xind < EnvXXXCfg.EnvWidth_c; xind++) {
        for (yind = 0; yind < EnvXXXCfg.EnvHeight_c; yind++) {
            EnvXXXCfg.Grid2D[xind][yind] = mapdata[xind + yind * EnvXXXCfg.EnvWidth_c];
        }
    }

    bNeedtoRecomputeStartHeuristics = true;
    bNeedtoRecomputeGoalHeuristics = true;

    return true;
}

void EnvironmentXXXLATTICE::PrintEnv_Config(FILE* fOut)
{
    // implement this if the planner needs to print out EnvXXXLAT. configuration
    throw SBPL_Exception("ERROR in EnvXXXLAT... function: PrintEnv_Config is undefined");
}

void EnvironmentXXXLATTICE::Set2DBlockSize(int BlockSize)
{
    blocksize = BlockSize;
}

void EnvironmentXXXLATTICE::Set2DBucketSize(int BucketSize)
{
    bucketsize = BucketSize;
}

void EnvironmentXXXLATTICE::PrintTimeStat(FILE* fOut)
{
#if TIME_DEBUG
    SBPL_FPRINTF(fOut, "time3_addallout = %f secs, time_gethash = %f secs, time_createhash = %f secs, "
                 "time_getsuccs = %f\n",
                 time3_addallout/(double)CLOCKS_PER_SEC, time_gethash/(double)CLOCKS_PER_SEC,
                 time_createhash/(double)CLOCKS_PER_SEC, time_getsuccs/(double)CLOCKS_PER_SEC);
#endif
}

bool EnvironmentXXXLATTICE::IsObstacle(int x, int y)
{
#if DEBUG
    SBPL_FPRINTF(fDeb, "Status of cell %d %d is queried. Its cost=%d\n", x,y,EnvXXXCfg.Grid2D[x][y]);
#endif

    return EnvXXXCfg.Grid2D[x][y] >= EnvXXXCfg.obsthresh;
}

void EnvironmentXXXLATTICE::GetEnvParms(
    int *size_x, int *size_y, int* num_thetas,
    double* startx, double* starty, double* starttheta,
    double* goalx, double* goaly, double* goaltheta,
    double* cellsize_m,
    double* nominalvel_mpersecs,
    double* timetoturn45degsinplace_secs,
    unsigned char* obsthresh,
    std::vector<SBPL_xxx_mprimitive>* mprimitiveV)
{
    *num_thetas = EnvXXXCfg.NumThetaDirs;
    GetEnvParms(
            size_x, size_y,
            startx, starty, starttheta,
            goalx, goaly, goaltheta,
            cellsize_m,
            nominalvel_mpersecs, timetoturn45degsinplace_secs,
            obsthresh,
            mprimitiveV);
}

void EnvironmentXXXLATTICE::GetEnvParms(
    int* size_x, int* size_y,
    double* startx, double* starty, double* starttheta,
    double* goalx, double* goaly, double* goaltheta,
    double* cellsize_m,
    double* nominalvel_mpersecs, double* timetoturn45degsinplace_secs,
    unsigned char* obsthresh,
    std::vector<SBPL_xxx_mprimitive>* mprimitiveV)
{
    *size_x = EnvXXXCfg.EnvWidth_c;
    *size_y = EnvXXXCfg.EnvHeight_c;

    *startx = DISCXY2CONT(EnvXXXCfg.StartX_c, EnvXXXCfg.cellsize_m);
    *starty = DISCXY2CONT(EnvXXXCfg.StartY_c, EnvXXXCfg.cellsize_m);
    *starttheta = DiscTheta2ContNew(EnvXXXCfg.StartTheta);
    *goalx = DISCXY2CONT(EnvXXXCfg.EndX_c, EnvXXXCfg.cellsize_m);
    *goaly = DISCXY2CONT(EnvXXXCfg.EndY_c, EnvXXXCfg.cellsize_m);
    *goaltheta = DiscTheta2ContNew(EnvXXXCfg.EndTheta);

    *cellsize_m = EnvXXXCfg.cellsize_m;
    *nominalvel_mpersecs = EnvXXXCfg.nominalvel_mpersecs;
    *timetoturn45degsinplace_secs = EnvXXXCfg.timetoturn45degsinplace_secs;

    *obsthresh = EnvXXXCfg.obsthresh;

    *mprimitiveV = EnvXXXCfg.mprimV;
}

bool EnvironmentXXXLATTICE::PoseContToDisc(
    double px, double py, double pth, int &ix, int &iy, int &ith) const
{
    ix = CONTXY2DISC(px, EnvXXXCfg.cellsize_m);
    iy = CONTXY2DISC(py, EnvXXXCfg.cellsize_m);
    ith = ContTheta2DiscNew(pth);
    return (pth >= -2 * PI_CONST) && (pth <= 2 * PI_CONST) &&
            (ix >= 0) && (ix < EnvXXXCfg.EnvWidth_c) &&
            (iy >= 0) && (iy < EnvXXXCfg.EnvHeight_c);
}

bool EnvironmentXXXLATTICE::PoseDiscToCont(
    int ix, int iy, int ith, double &px, double &py, double &pth) const
{
    px = DISCXY2CONT(ix, EnvXXXCfg.cellsize_m);
    py = DISCXY2CONT(iy, EnvXXXCfg.cellsize_m);
    pth = normalizeAngle(DiscTheta2ContNew(ith));
    return (ith >= 0) && (ith < EnvXXXCfg.NumThetaDirs) &&
            (ix >= 0) && (ix < EnvXXXCfg.EnvWidth_c) &&
            (iy >= 0) && (iy < EnvXXXCfg.EnvHeight_c);
}

unsigned char EnvironmentXXXLATTICE::GetMapCost(int x, int y)
{
    return EnvXXXCfg.Grid2D[x][y];
}

bool EnvironmentXXXLATTICE::SetEnvParameter(
    const char* parameter,
    int value)
{
    if (EnvXXXLAT.bInitialized) {
        SBPL_ERROR("ERROR: all parameters must be set before initialization of the environment\n");
        return false;
    }

    SBPL_PRINTF("setting parameter %s to %d\n", parameter, value);

    if (strcmp(parameter, "cost_inscribed_thresh") == 0) {
        if (value < 0 || value > 255) {
            SBPL_ERROR("ERROR: invalid value %d for parameter %s\n", value, parameter);
            return false;
        }
        EnvXXXCfg.cost_inscribed_thresh = (unsigned char)value;
    }
    else if (strcmp(parameter, "cost_possibly_circumscribed_thresh") == 0) {
        if (value < 0 || value > 255) {
            SBPL_ERROR("ERROR: invalid value %d for parameter %s\n", value, parameter);
            return false;
        }
        EnvXXXCfg.cost_possibly_circumscribed_thresh = value;
    }
    else if (strcmp(parameter, "cost_obsthresh") == 0) {
        if (value < 0 || value > 255) {
            SBPL_ERROR("ERROR: invalid value %d for parameter %s\n", value, parameter);
            return false;
        }
        EnvXXXCfg.obsthresh = (unsigned char)value;
    }
    else {
        SBPL_ERROR("ERROR: invalid parameter %s\n", parameter);
        return false;
    }

    return true;
}

int EnvironmentXXXLATTICE::GetEnvParameter(const char* parameter)
{
    if (strcmp(parameter, "cost_inscribed_thresh") == 0) {
        return (int)EnvXXXCfg.cost_inscribed_thresh;
    }
    else if (strcmp(parameter, "cost_possibly_circumscribed_thresh") == 0) {
        return (int)EnvXXXCfg.cost_possibly_circumscribed_thresh;
    }
    else if (strcmp(parameter, "cost_obsthresh") == 0) {
        return (int)EnvXXXCfg.obsthresh;
    }
    else {
        std::stringstream ss;
        ss << "ERROR: invalid parameter " << parameter;
        throw SBPL_Exception(ss.str());
    }
}

EnvironmentXXXLAT::~EnvironmentXXXLAT()
{
    SBPL_PRINTF("destroying XYTHETALAT\n");

    // delete the states themselves first
    for (int i = 0; i < (int)StateID2CoordTable.size(); i++) {
        delete StateID2CoordTable.at(i);
        StateID2CoordTable.at(i) = NULL;
    }
    StateID2CoordTable.clear();
    stateToPathHistory.clear();

    // delete hashtable
    if (Coord2StateIDHashTable != NULL) {
        delete[] Coord2StateIDHashTable;
        Coord2StateIDHashTable = NULL;
    }
    if (Coord2StateIDHashTable_lookup != NULL) {
        delete[] Coord2StateIDHashTable_lookup;
        Coord2StateIDHashTable_lookup = NULL;
    }
}

void EnvironmentXXXLAT::GetCoordFromState(
    int stateID, int& x, int& y, int& theta) const
{
    EnvXXXLATHashEntry_t* HashEntry = StateID2CoordTable[stateID];
    x = HashEntry->X;
    y = HashEntry->Y;
    theta = HashEntry->Theta;
}

int EnvironmentXXXLAT::GetStateFromCoord(int x, int y, int theta)
{
    EnvXXXLATHashEntry_t* OutHashEntry;
    if ((OutHashEntry = (this->*GetHashEntry)(x, y, theta)) == NULL) {
        // have to create a new entry
        OutHashEntry = (this->*CreateNewHashEntry)(x, y, theta);
    }
    return OutHashEntry->stateID;
}

void EnvironmentXXXLAT::GetActionsFromStateIDPath(
    std::vector<int>* stateIDPath,
    std::vector<EnvXXXLATAction_t>* action_list)
{
    std::vector<EnvXXXLATAction_t*> actionV;
    std::vector<int> CostV;
    std::vector<int> SuccIDV;
    int targetx_c, targety_c, targettheta_c;
    int sourcex_c, sourcey_c, sourcetheta_c;

    SBPL_PRINTF("checks=%ld\n", checks);

    action_list->clear();

    for (int pind = 0; pind < (int)(stateIDPath->size()) - 1; pind++) {
        int sourceID = stateIDPath->at(pind);
        int targetID = stateIDPath->at(pind + 1);

        // get successors and pick the target via the cheapest action
        SuccIDV.clear();
        CostV.clear();
        actionV.clear();
        GetSuccs(sourceID, &SuccIDV, &CostV, &actionV);

        int bestcost = INFINITECOST;
        int bestsind = -1;

        for (int sind = 0; sind < (int)SuccIDV.size(); sind++) {
            if (SuccIDV[sind] == targetID && CostV[sind] <= bestcost) {
                bestcost = CostV[sind];
                bestsind = sind;
            }
        }
        if (bestsind == -1) {
            SBPL_ERROR("ERROR: successor not found for transition");
            GetCoordFromState(sourceID, sourcex_c, sourcey_c, sourcetheta_c);
            GetCoordFromState(targetID, targetx_c, targety_c, targettheta_c);
            SBPL_PRINTF("%d %d %d -> %d %d %d\n", sourcex_c, sourcey_c, sourcetheta_c, targetx_c, targety_c, targettheta_c);
            throw SBPL_Exception("ERROR: successor not found for transition");
        }

#if DEBUG
        SBPL_FPRINTF(fDeb, "Start: %.3f %.3f %.3f Target: %.3f %.3f %.3f Prim ID, Start Theta: %d %d\n",
            sourcex_c, sourcey_c, sourcetheta_c,
            targetx_c, targety_c, targettheta_c,
            actionV[bestsind]->aind, actionV[bestsind]->starttheta);
#endif

        action_list->push_back(*(actionV[bestsind]));
    }
}

void EnvironmentXXXLAT::ConvertStateIDPathintoXYThetaPath(
    std::vector<int>* stateIDPath,
    std::vector<sbpl_xy_theta_pt_t>* xythetaPath,
    std::vector<sbpl_xy_theta_pt_t>* trailerPath)
{
    std::vector<EnvXXXLATAction_t*> actionV;
    std::vector<int> CostV;
    std::vector<int> SuccIDV;
    int targetx_c, targety_c, targettheta_c;
    int sourcex_c, sourcey_c, sourcetheta_c;

    SBPL_PRINTF("checks=%ld\n", checks);

    xythetaPath->clear();
    trailerPath->clear();

    // std::ofstream logFile("/home/opencav/Downloads/sbpl_trailerpose.log", std::ios::app);
    // if (!logFile.is_open()) {
    //     throw SBPL_Exception("ERROR: Could not open log for writing");
    // }

#if DEBUG
    SBPL_FPRINTF(fDeb, "converting stateid path into coordinates:\n");
#endif

    for (int pind = 0; pind < (int)(stateIDPath->size()) - 1; pind++) {
        int sourceID = stateIDPath->at(pind);
        int targetID = stateIDPath->at(pind + 1);

        // auto sourceHistoryIt = stateToPathHistory.find(sourceID);
        // if (sourceHistoryIt == stateToPathHistory.end()) {
        //     SBPL_ERROR("Path history not fund for state %d", sourceID);
        //     throw SBPL_Exception("ERROR: path history not found");
        // }
        // const PathHistory& currentHistory = sourceHistoryIt->second;
        // TrailerState currentTrailer = currentHistory.currentTrailer;
        auto targetHistoryIt = stateToPathHistory.find(targetID);
        if (targetHistoryIt == stateToPathHistory.end()) {
            SBPL_ERROR("Path history not fund for state %d", targetID);
            throw SBPL_Exception("ERROR: path history not found");
        }
        const PathHistory& currentHistory = targetHistoryIt->second;
        TrailerState currentTrailer = currentHistory.currentTrailer;

#if DEBUG
        GetCoordFromState(sourceID, sourcex_c, sourcey_c, sourcetheta_c);
#endif

        // get successors and pick the target via the cheapest action
        SuccIDV.clear();
        CostV.clear();
        actionV.clear();
        // GetSuccs(sourceID, &SuccIDV, &CostV, &actionV);
        GetSuccsOfBestPath(sourceID, &SuccIDV, &CostV, &actionV);

        int bestcost = INFINITECOST;
        int bestsind = -1;

#if DEBUG
        GetCoordFromState(sourceID, sourcex_c, sourcey_c, sourcetheta_c);
        GetCoordFromState(targetID, targetx_c, targety_c, targettheta_c);
        SBPL_FPRINTF(fDeb, "looking for %d %d %d -> %d %d %d (numofsuccs=%d)\n", sourcex_c, sourcey_c, sourcetheta_c, targetx_c, targety_c, targettheta_c, (int)SuccIDV.size());
#endif

        for (int sind = 0; sind < (int)SuccIDV.size(); sind++) {
#if DEBUG
            int x_c, y_c, theta_c;
            GetCoordFromState(SuccIDV[sind], x_c, y_c, theta_c);
            SBPL_FPRINTF(fDeb, "succ: %d %d %d\n", x_c, y_c, theta_c);
#endif
            if (SuccIDV[sind] == targetID && CostV[sind] <= bestcost) {
                bestcost = CostV[sind];
                bestsind = sind;
            }
        }
        if (bestsind == -1) {
            SBPL_ERROR("ERROR: successor not found for transition");
            GetCoordFromState(sourceID, sourcex_c, sourcey_c, sourcetheta_c);
            GetCoordFromState(targetID, targetx_c, targety_c, targettheta_c);
            SBPL_PRINTF("%d %d %d -> %d %d %d\n", sourcex_c, sourcey_c, sourcetheta_c, targetx_c, targety_c, targettheta_c);
            SBPL_PRINTF("%.3f %.3f %.3f -> %.3f %.3f %.3f\n", DISCXY2CONT(sourcex_c, EnvXXXCfg.cellsize_m), DISCXY2CONT(sourcey_c, EnvXXXCfg.cellsize_m), DiscTheta2ContNew(sourcetheta_c), DISCXY2CONT(targetx_c, EnvXXXCfg.cellsize_m), DISCXY2CONT(targety_c, EnvXXXCfg.cellsize_m), DiscTheta2ContNew(targettheta_c));
            // throw SBPL_WARN("ERROR: successor not found for transition");
        }

        // now push in the actual path
        GetCoordFromState(sourceID, sourcex_c, sourcey_c, sourcetheta_c);
        double sourcex, sourcey, sourcetheta;
        sourcex = DISCXY2CONT(sourcex_c, EnvXXXCfg.cellsize_m);
        sourcey = DISCXY2CONT(sourcey_c, EnvXXXCfg.cellsize_m);
        sourcetheta = DiscTheta2ContNew(sourcetheta_c);

        // double currentX = sourcex;
        // double currentY = sourcey;
        // double currentTheta = sourcetheta;

        // std::vector<PathState> intermediatePath = currentHistory.positions;

        // TODO - when there are no motion primitives we should still print source state
        for (int ipind = 0; ipind < ((int)actionV[bestsind]->intermptV.size()) - 1; ipind++) {
            // translate appropriately
            sbpl_xy_theta_pt_t intermpt = actionV[bestsind]->intermptV[ipind];
            intermpt.x += sourcex;
            intermpt.y += sourcey;

#if DEBUG
            int nx = CONTXY2DISC(intermpt.x, EnvXXXCfg.cellsize_m);
            int ny = CONTXY2DISC(intermpt.y, EnvXXXCfg.cellsize_m);
            int ntheta;
            ntheta = ContTheta2DiscNew(intermpt.theta);

            SBPL_FPRINTF(fDeb, "%.3f %.3f %.3f (%d %d %d cost=%d) ", intermpt.x, intermpt.y, intermpt.theta, nx, ny, ntheta,  EnvXXXCfg.Grid2D[nx][ny]);

            if (ipind == 0) {
                SBPL_FPRINTF(fDeb, "first (heur=%d)\n", GetStartHeuristic(sourceID));
            }
            else {
                SBPL_FPRINTF(fDeb, "\n");
            }
#endif
            // store
            xythetaPath->push_back(intermpt);

            // logFile << "Step " << ipind << " of primitive " << bestsind << "\n";
            // logFile << "    Source State: x=" << currentX << ", y=" << currentY << ", theta=" << currentTheta << "\n";
            // logFile << "    Target State: x=" << intermpt.x << ", y=" << intermpt.y << ", theta=" << intermpt.theta << "\n";
            // currentX = intermpt.x;
            // currentY = intermpt.y;
            // currentTheta = intermpt.theta;
        }
        // PathState intermediateState = {
        //     currentX, currentY, currentTheta,
        //     currentHistory.positions.back().time + actionV[bestsind]->time
        // };
        // intermediatePath.push_back(intermediateState);
        // TrailerState intermediateTrailer;
        // if (calculateTrailerFromPath(intermediatePath, intermediateTrailer)) {
        //     // logFile << sourceID << "- ROBOT values -- x: (" << intermediateState.x << " " << intermpt.x << "), y: (" << intermediateState.y << " " << intermpt.y << "), theta: (" << intermediateState.theta << " " << intermpt.theta << ")" << "\n";
        //     logFile << sourceID << "- TRAILER values -- x: (" << currentTrailer.x << " " << intermediateTrailer.x << "), y: (" << currentTrailer.y << " " << intermediateTrailer.y << "), theta1: (" << currentTrailer.theta1 << " " << intermediateTrailer.theta1 << "), theta2: (" << currentTrailer.theta2 << " " << intermediateTrailer.theta2 << ")" << "\n";
        //     logFile << "\n";
        //     sbpl_xy_theta_pt_t trailer_pt;
        //     trailer_pt.x = intermediateTrailer.x;
        //     trailer_pt.y = intermediateTrailer.y;
        //     trailer_pt.theta = intermediateTrailer.theta2;
        //     trailerPath->push_back(trailer_pt);
        //     currentTrailer = intermediateTrailer;
        // }
        sbpl_xy_theta_pt_t trailer_pt;
        trailer_pt.x = currentTrailer.x;
        trailer_pt.y = currentTrailer.y;
        trailer_pt.theta = currentTrailer.theta2;
        trailerPath->push_back(trailer_pt);

    }
}

// returns the stateid if success, and -1 otherwise
int EnvironmentXXXLAT::SetGoal(double x_m, double y_m, double theta_rad)
{
    int x = CONTXY2DISC(x_m, EnvXXXCfg.cellsize_m);
    int y = CONTXY2DISC(y_m, EnvXXXCfg.cellsize_m);
    int theta = ContTheta2DiscNew(theta_rad);

    SBPL_PRINTF("env: setting goal to %.3f %.3f %.3f (%d %d %d)\n", x_m, y_m, theta_rad, x, y, theta);

    if (!IsWithinMapCell(x, y)) {
        SBPL_ERROR("ERROR: trying to set a goal cell %d %d that is outside of map\n", x, y);
        return -1;
    }

    if (!IsValidConfiguration(x, y, theta)) {
        SBPL_PRINTF("WARNING: goal configuration is invalid\n");
    }

    EnvXXXLATHashEntry_t* OutHashEntry;
    if ((OutHashEntry = (this->*GetHashEntry)(x, y, theta)) == NULL) {
        // have to create a new entry
        OutHashEntry = (this->*CreateNewHashEntry)(x, y, theta);
    }

    // need to recompute start heuristics?
    if (EnvXXXLAT.goalstateid != OutHashEntry->stateID) {
        // because termination condition may not plan all the way to the new goal
        bNeedtoRecomputeStartHeuristics = true;

        // because goal heuristics change
        bNeedtoRecomputeGoalHeuristics = true;
    }

    EnvXXXLAT.goalstateid = OutHashEntry->stateID;

    EnvXXXCfg.EndX_c = x;
    EnvXXXCfg.EndY_c = y;
    EnvXXXCfg.EndTheta = theta;

    return EnvXXXLAT.goalstateid;
}

// returns the stateid if success, and -1 otherwise
int EnvironmentXXXLAT::SetStart(double x_m, double y_m, double theta_rad)
{
    int x = CONTXY2DISC(x_m, EnvXXXCfg.cellsize_m);
    int y = CONTXY2DISC(y_m, EnvXXXCfg.cellsize_m);
    int theta = ContTheta2DiscNew(theta_rad);

    if (!IsWithinMapCell(x, y)) {
        SBPL_ERROR("ERROR: trying to set a start cell %d %d that is outside of map\n", x, y);
        return -1;
    }

    SBPL_PRINTF("env: setting start to %.3f %.3f %.3f (%d %d %d)\n", x_m, y_m, theta_rad, x, y, theta);

    if (!IsValidConfiguration(x, y, theta)) {
        SBPL_PRINTF("WARNING: start configuration %d %d %d is invalid\n", x, y, theta);
    }

    EnvXXXLATHashEntry_t* OutHashEntry;
    if ((OutHashEntry = (this->*GetHashEntry)(x, y, theta)) == NULL) {
        // have to create a new entry
        OutHashEntry = (this->*CreateNewHashEntry)(x, y, theta);
    }

    // need to recompute start heuristics?
    if (EnvXXXLAT.startstateid != OutHashEntry->stateID) {
        bNeedtoRecomputeStartHeuristics = true;
        // because termination condition can be not all states TODO - make it dependent on term. condition
        bNeedtoRecomputeGoalHeuristics = true;
    }

    // set start
    EnvXXXLAT.startstateid = OutHashEntry->stateID;
    EnvXXXCfg.StartX_c = x;
    EnvXXXCfg.StartY_c = y;
    EnvXXXCfg.StartTheta = theta;

    return EnvXXXLAT.startstateid;
}

// returns the stateid if success, and -1 otherwise
void EnvironmentXXXLAT::SetTrailerStart(double theta1_rad, double theta2_rad)
{
    EnvXXXCfg.StartTheta1 = theta1_rad;
    EnvXXXCfg.StartTheta2 = theta2_rad;
}

void EnvironmentXXXLAT::PrintState(
    int stateID,
    bool bVerbose,
    FILE* fOut)
{
#if DEBUG
    if (stateID >= (int)StateID2CoordTable.size()) {
        SBPL_ERROR("ERROR in EnvXXXLAT... function: stateID illegal (2)\n");
        throw SBPL_Exception();
    }
#endif

    if (fOut == NULL) {
        fOut = stdout;
    }

    EnvXXXLATHashEntry_t* HashEntry = StateID2CoordTable[stateID];

    if (stateID == EnvXXXLAT.goalstateid && bVerbose) {
        SBPL_FPRINTF(fOut, "the state is a goal state\n");
    }

    if (bVerbose) {
        SBPL_FPRINTF(fOut, "X=%d Y=%d Theta=%d\n", HashEntry->X, HashEntry->Y, HashEntry->Theta);
    }
    else
    {
        SBPL_FPRINTF(fOut, "%.3f %.3f %.3f\n", DISCXY2CONT(HashEntry->X, EnvXXXCfg.cellsize_m), DISCXY2CONT(HashEntry->Y, EnvXXXCfg.cellsize_m), DiscTheta2ContNew(HashEntry->Theta));
    }
}

EnvXXXLATHashEntry_t* EnvironmentXXXLAT::GetHashEntry_lookup(
    int X, int Y, int Theta)
{
    if (X < 0 || X >= EnvXXXCfg.EnvWidth_c ||
        Y < 0 || Y >= EnvXXXCfg.EnvHeight_c ||
        Theta < 0 || Theta >= EnvXXXCfg.NumThetaDirs)
    {
        return NULL;
    }
    int index = XYTHETA2INDEX(X,Y,Theta);
    return Coord2StateIDHashTable_lookup[index];
}

EnvXXXLATHashEntry_t*
EnvironmentXXXLAT::GetHashEntry_hash(int X, int Y, int Theta)
{
#if TIME_DEBUG
    clock_t currenttime = clock();
#endif

    int binid = GETHASHBIN(X, Y, Theta);

#if DEBUG
    if ((int)Coord2StateIDHashTable[binid].size() > 5) {
        SBPL_FPRINTF(fDeb, "WARNING: Hash table has a bin %d (X=%d Y=%d) of size %d\n", binid, X, Y, (int)Coord2StateIDHashTable[binid].size());

        PrintHashTableHist(fDeb);
    }
#endif

    // iterate over the states in the bin and select the perfect match
    std::vector<EnvXXXLATHashEntry_t*>* binV = &Coord2StateIDHashTable[binid];
    for (int ind = 0; ind < (int)binV->size(); ind++) {
        EnvXXXLATHashEntry_t* hashentry = binV->at(ind);
        if (hashentry->X == X && hashentry->Y == Y && hashentry->Theta == Theta) {
#if TIME_DEBUG
            time_gethash += clock()-currenttime;
#endif
            return hashentry;
        }
    }

#if TIME_DEBUG
    time_gethash += clock()-currenttime;
#endif

    return NULL;
}

EnvXXXLATHashEntry_t*
EnvironmentXXXLAT::CreateNewHashEntry_lookup(int X, int Y, int Theta)
{
    int i;

#if TIME_DEBUG
    clock_t currenttime = clock();
#endif

    EnvXXXLATHashEntry_t* HashEntry = new EnvXXXLATHashEntry_t;

    HashEntry->X = X;
    HashEntry->Y = Y;
    HashEntry->Theta = Theta;
    HashEntry->iteration = 0;

    HashEntry->stateID = StateID2CoordTable.size();

    // insert into the tables
    StateID2CoordTable.push_back(HashEntry);

    int index = XYTHETA2INDEX(X,Y,Theta);

#if DEBUG
    if (Coord2StateIDHashTable_lookup[index] != NULL) {
        throw SBPL_Exception("ERROR: creating hash entry for non-NULL hashentry");
    }
#endif

    Coord2StateIDHashTable_lookup[index] = HashEntry;

    // insert into and initialize the mappings
    int* entry = new int[NUMOFINDICES_STATEID2IND];
    StateID2IndexMapping.push_back(entry);
    for (i = 0; i < NUMOFINDICES_STATEID2IND; i++) {
        StateID2IndexMapping[HashEntry->stateID][i] = -1;
    }

    if (HashEntry->stateID != (int)StateID2IndexMapping.size() - 1) {
        throw SBPL_Exception("ERROR in Env... function: last state has incorrect stateID");
    }

#if TIME_DEBUG
    time_createhash += clock()-currenttime;
#endif

    return HashEntry;
}

EnvXXXLATHashEntry_t*
EnvironmentXXXLAT::CreateNewHashEntry_hash(int X, int Y, int Theta)
{
    int i;

#if TIME_DEBUG
    clock_t currenttime = clock();
#endif

    EnvXXXLATHashEntry_t* HashEntry = new EnvXXXLATHashEntry_t;

    HashEntry->X = X;
    HashEntry->Y = Y;
    HashEntry->Theta = Theta;
    HashEntry->iteration = 0;

    HashEntry->stateID = StateID2CoordTable.size();

    // insert into the tables
    StateID2CoordTable.push_back(HashEntry);

    // get the hash table bin
    i = GETHASHBIN(HashEntry->X, HashEntry->Y, HashEntry->Theta);

    // insert the entry into the bin
    Coord2StateIDHashTable[i].push_back(HashEntry);

    // insert into and initialize the mappings
    int* entry = new int[NUMOFINDICES_STATEID2IND];
    StateID2IndexMapping.push_back(entry);
    for (i = 0; i < NUMOFINDICES_STATEID2IND; i++) {
        StateID2IndexMapping[HashEntry->stateID][i] = -1;
    }

    if (HashEntry->stateID != (int)StateID2IndexMapping.size() - 1) {
        throw SBPL_Exception("ERROR in Env... function: last state has incorrect stateID");
    }

#if TIME_DEBUG
    time_createhash += clock() - currenttime;
#endif

    return HashEntry;
}

double EnvironmentXXXLAT::calculateContAngleDiff(double startangle, double endangle){
    double diff = endangle - startangle;
    if (diff <= -PI_CONST){
        diff += 2*PI_CONST;
    } else if (diff > PI_CONST){
        diff -= 2*PI_CONST;
    }
    return diff;
}

bool EnvironmentXXXLATTICE::calculateTrailerFromPath(
    const std::vector<PathState>& path,
    TrailerState& finalTrailer) 
{
    if (path.size() < 2) return false;

    // TrailerState trailer(
    //     path[0].x - R0*cos(path[0].theta) - F1*cos(path[0].theta) - (F2/2)*cos(path[0].theta),
    //     path[0].y - R0*sin(path[0].theta) - F1*sin(path[0].theta) - (F2/2)*sin(path[0].theta),
    //     path[0].theta,
    //     path[0].theta
    // );
    TrailerState trailer(
        path[0].x - R0*cos(path[0].theta) - F1*cos(EnvXXXCfg.StartTheta1) - (F2/2)*cos(EnvXXXCfg.StartTheta2),
        path[0].y - R0*sin(path[0].theta) - F1*sin(EnvXXXCfg.StartTheta1) - (F2/2)*sin(EnvXXXCfg.StartTheta2),
        EnvXXXCfg.StartTheta1,
        EnvXXXCfg.StartTheta2
    );

    for (int i = 1; i < path.size(); i++) {
        const PathState& prev = path[i-1];
        const PathState& curr = path[i];
        double deltaTime = curr.time - prev.time;
        
        if (!calculateTrailerTransition(
            prev.x, prev.y, prev.theta,
            curr.x, curr.y, curr.theta,
            deltaTime,
            trailer,
            trailer)) //update trailer in place
        {
            return false;
        }
    }
    finalTrailer = trailer;
    return true;
}

bool EnvironmentXXXLATTICE::calculateTrailerTransition(double startX, double startY, double startTheta, double endX, double endY, double endTheta,
                                                    double time, const TrailerState& startTrailer, TrailerState& endTrailer) {
    double t = time;
    const int steps = 10;
    double dt = t / steps;
    double v = EnvXXXCfg.nominalvel_mpersecs;

    double dx = endX - startX;
    double dy = endY - startY;
    double dtheta = calculateContAngleDiff(startTheta, endTheta);
    double w = dtheta / t;

    double theta = startTheta;
    double theta1 = startTrailer.theta1;
    double theta2 = startTrailer.theta2;

    for (int i = 0; i < steps; i++) {
        double phi1 = calculateContAngleDiff(theta1, theta);
        double phi2 = calculateContAngleDiff(theta2, theta1);

        double v1 = v*cos(phi1) + R0*w*sin(phi1);
        double dtheta1 = (v*sin(phi1) - R0*w*cos(phi1))/F1;
        double dtheta2 = v1*sin(phi2)/F2;

        theta += w*dt;
        theta1 += dtheta1*dt;
        theta2 += dtheta2*dt;

        // theta = normalizeAngle(theta);
        // theta1 = normalizeAngle(theta1);
        // theta2 = normalizeAngle(theta2);
    }

    endTrailer.theta1 = theta1;
    endTrailer.theta2 = theta2;
    endTrailer.x = endX - R0*cos(endTheta) - F1*cos(theta1) - (F2/2)*cos(theta2);
    endTrailer.y = endY - R0*sin(endTheta) - F1*sin(theta1) - (F2/2)*sin(theta2);

    // int disc_endX = CONTXY2DISC(endTrailer.x, EnvXXXCfg.cellsize_m);
    // int disc_endY = CONTXY2DISC(endTrailer.y, EnvXXXCfg.cellsize_m);
    // if (!IsValidCell(disc_endX, disc_endY)) {
    //     return false;
    // }
    // if (EnvXXXCfg.Grid2D[disc_endX][disc_endY] >= EnvXXXCfg.cost_inscribed_thresh)
    // {
    //     return false;
    // }
    // return true;
    return IsValidTrailerConfiguration(endTrailer.x, endTrailer.y, endTrailer.theta2);
}

bool EnvironmentXXXLATTICE::IsValidTrailerConfiguration(double x, double y, double theta) {
    int disc_x = CONTXY2DISC(x, EnvXXXCfg.cellsize_m);
    int disc_y = CONTXY2DISC(y, EnvXXXCfg.cellsize_m);
    if (!IsValidCell(disc_x, disc_y)) {
        return false;
    }
    if (EnvXXXCfg.Grid2D[disc_x][disc_y] >= EnvXXXCfg.cost_inscribed_thresh) {
        return false;
    }
    
    std::vector<sbpl_2Dpt_t> transformedPolygon;
    transformedPolygon.resize(EnvXXXCfg.TrailerPolygon.size());

    for (int i = 0; i < EnvXXXCfg.TrailerPolygon.size(); i++) {
        transformedPolygon[i].x = x + EnvXXXCfg.TrailerPolygon[i].x * cos(theta) - EnvXXXCfg.TrailerPolygon[i].y * sin(theta);
        transformedPolygon[i].y = y + EnvXXXCfg.TrailerPolygon[i].x * sin(theta) + EnvXXXCfg.TrailerPolygon[i].y * cos(theta);
    }

    for (const auto& pt : transformedPolygon) {
        int poly_x = CONTXY2DISC(pt.x, EnvXXXCfg.cellsize_m);
        int poly_y = CONTXY2DISC(pt.y, EnvXXXCfg.cellsize_m);
        if (!IsValidCell(poly_x, poly_y)) {
            return false;
        }
        if (EnvXXXCfg.Grid2D[poly_x][poly_y] >= EnvXXXCfg.cost_inscribed_thresh) {
            return false;
        }
    }
    return true;
}

void EnvironmentXXXLAT::GetSuccs(
    int SourceStateID,
    std::vector<int>* SuccIDV,
    std::vector<int>* CostV,
    std::vector<EnvXXXLATAction_t*>* actionV)
{
    int aind;

#if TIME_DEBUG
    clock_t currenttime = clock();
#endif

    // clear the successor array
    SuccIDV->clear();
    CostV->clear();
    SuccIDV->reserve(EnvXXXCfg.actionwidth);
    CostV->reserve(EnvXXXCfg.actionwidth);
    if (actionV != NULL) {
        actionV->clear();
        actionV->reserve(EnvXXXCfg.actionwidth);
    }

    // goal state should be absorbing
    if (SourceStateID == EnvXXXLAT.goalstateid) return;

    // get X, Y for the state
    EnvXXXLATHashEntry_t* HashEntry = StateID2CoordTable[SourceStateID];

    double currentX, currentY, currentTheta;
    currentX = DISCXY2CONT(HashEntry->X, EnvXXXCfg.cellsize_m);
    currentY = DISCXY2CONT(HashEntry->Y, EnvXXXCfg.cellsize_m);
    currentTheta = DiscTheta2ContNew(HashEntry->Theta);

    PathHistory currentHistory;
    TrailerState currentTrailer;
    // If we are at the start then set the trailer pose to theta1 = 0, theta2 = 0
    if (SourceStateID == EnvXXXLAT.startstateid) {
        PathState startState = {
            currentX,
            currentY,
            currentTheta,
            0.0
        };
        currentHistory.positions.push_back(startState);
        currentHistory.currentTrailer = TrailerState(
            currentX - R0*cos(currentTheta) - F1*cos(EnvXXXCfg.StartTheta1) - (F2/2)*cos(EnvXXXCfg.StartTheta2),
            currentY - R0*sin(currentTheta) - F1*sin(EnvXXXCfg.StartTheta1) - (F2/2)*sin(EnvXXXCfg.StartTheta2),
            EnvXXXCfg.StartTheta1,
            EnvXXXCfg.StartTheta2
        );
        stateToPathHistory[SourceStateID] = currentHistory;
        // trailerStates[SourceStateID] = currentTrailer;
    } else {
        auto it = stateToPathHistory.find(SourceStateID);
        if (it == stateToPathHistory.end()) {
            SBPL_ERROR("Path history not found for StateID %d", SourceStateID);
            return;
        }
        currentHistory = it->second;
    }

    // iterate through actions
    for (aind = 0; aind < EnvXXXCfg.actionwidth; aind++) {
        EnvXXXLATAction_t* nav3daction = &EnvXXXCfg.ActionsV[(unsigned int)HashEntry->Theta][aind];
        int newX = HashEntry->X + nav3daction->dX;
        int newY = HashEntry->Y + nav3daction->dY;
        int newTheta = normalizeDiscAngle(nav3daction->endtheta);
        double time = nav3daction->time;

        // skip the invalid cells
        if (!IsValidCell(newX, newY)) {
            continue;
        }

        // get cost
        int cost = GetActionCost(HashEntry->X, HashEntry->Y, HashEntry->Theta, nav3daction);
        if (cost >= INFINITECOST) {
            continue;
        }

        double cont_newX, cont_newY, cont_newTheta;
        cont_newX = DISCXY2CONT(newX, EnvXXXCfg.cellsize_m);
        cont_newY = DISCXY2CONT(newY, EnvXXXCfg.cellsize_m);
        cont_newTheta = DiscTheta2ContNew(newTheta);

        PathHistory newHistory = currentHistory;
        PathState newState = {
            cont_newX,
            cont_newY,
            cont_newTheta,
            currentHistory.positions.back().time + time
        };
        newHistory.positions.push_back(newState);

        TrailerState newTrailer;
        if (!calculateTrailerFromPath(newHistory.positions, newTrailer)) {
            continue;
        }
        newHistory.currentTrailer = newTrailer;
        
        EnvXXXLATHashEntry_t* OutHashEntry;
        if ((OutHashEntry = (this->*GetHashEntry)(newX, newY, newTheta)) == NULL) {
            // have to create a new entry
            OutHashEntry = (this->*CreateNewHashEntry)(newX, newY, newTheta);
        }

        stateToPathHistory[OutHashEntry->stateID] = std::move(newHistory);
        // SetStateParent(OutHashEntry->stateID, SourceStateID);
        // StateTransitionKey key{SourceStateID, OutHashEntry->stateID};
        // trailerTransitions[key] = newTrailer;
        // trailerStates[OutHashEntry->stateID] = newTrailer;

        SuccIDV->push_back(OutHashEntry->stateID);
        CostV->push_back(cost);
        if (actionV != NULL) {
            actionV->push_back(nav3daction);
        }
    }

#if TIME_DEBUG
    time_getsuccs += clock()-currenttime;
#endif
}

void EnvironmentXXXLAT::GetSuccsOfBestPath(
    int SourceStateID,
    std::vector<int>* SuccIDV,
    std::vector<int>* CostV,
    std::vector<EnvXXXLATAction_t*>* actionV)
{
    int aind;

    // clear the successor array
    SuccIDV->clear();
    CostV->clear();
    SuccIDV->reserve(EnvXXXCfg.actionwidth);
    CostV->reserve(EnvXXXCfg.actionwidth);
    if (actionV != NULL) {
        actionV->clear();
        actionV->reserve(EnvXXXCfg.actionwidth);
    }

    // goal state should be absorbing
    if (SourceStateID == EnvXXXLAT.goalstateid) return;

    // get X, Y for the state
    EnvXXXLATHashEntry_t* HashEntry = StateID2CoordTable[SourceStateID];

    // iterate through actions
    for (aind = 0; aind < EnvXXXCfg.actionwidth; aind++) {
        EnvXXXLATAction_t* nav3daction = &EnvXXXCfg.ActionsV[(unsigned int)HashEntry->Theta][aind];
        int newX = HashEntry->X + nav3daction->dX;
        int newY = HashEntry->Y + nav3daction->dY;
        int newTheta = normalizeDiscAngle(nav3daction->endtheta);

        // SBPL_DEBUG("Checking successor: current(%d,%d) -> new(%d,%d)", HashEntry->X, HashEntry->Y, newX, newY);

        // skip the invalid cells
        if (!IsValidCell(newX, newY)) {
            SBPL_DEBUG("Invalid cell at (%d,%d)", newX, newY);
            continue;
        }

        // get cost
        int cost = GetActionCost(HashEntry->X, HashEntry->Y, HashEntry->Theta, nav3daction);
        if (cost >= INFINITECOST) {
            SBPL_DEBUG("Infinite cost for action %d", aind);
            continue;
        }
        
        EnvXXXLATHashEntry_t* OutHashEntry;
        if ((OutHashEntry = (this->*GetHashEntry)(newX, newY, newTheta)) == NULL) {
            // have to create a new entry
            OutHashEntry = (this->*CreateNewHashEntry)(newX, newY, newTheta);
        }

        SuccIDV->push_back(OutHashEntry->stateID);
        CostV->push_back(cost);
        if (actionV != NULL) {
            actionV->push_back(nav3daction);
        }
    }
}

void EnvironmentXXXLAT::GetPreds(
    int TargetStateID,
    std::vector<int>* PredIDV,
    std::vector<int>* CostV)
{
    //TODO- to support tolerance, need:
    // a) generate preds for goal state based on all possible goal state variable settings,
    // b) change goal check condition in gethashentry c) change
    //    getpredsofchangedcells and getsuccsofchangedcells functions

    int aind;

#if TIME_DEBUG
    clock_t currenttime = clock();
#endif

    // get X, Y for the state
    EnvXXXLATHashEntry_t* HashEntry = StateID2CoordTable[TargetStateID];

    // clear the successor array
    PredIDV->clear();
    CostV->clear();
    PredIDV->reserve(EnvXXXCfg.PredActionsV[(unsigned int)HashEntry->Theta].size());
    CostV->reserve(EnvXXXCfg.PredActionsV[(unsigned int)HashEntry->Theta].size());

    // iterate through actions
    std::vector<EnvXXXLATAction_t*>* actionsV = &EnvXXXCfg.PredActionsV[(unsigned int)HashEntry->Theta];
    for (aind = 0; aind < (int)EnvXXXCfg.PredActionsV[(unsigned int)HashEntry->Theta].size(); aind++) {

        EnvXXXLATAction_t* nav3daction = actionsV->at(aind);

        int predX = HashEntry->X - nav3daction->dX;
        int predY = HashEntry->Y - nav3daction->dY;
        int predTheta = nav3daction->starttheta;

        // skip the invalid cells
        if (!IsValidCell(predX, predY)) {
            continue;
        }

        // get cost
        int cost = GetActionCost(predX, predY, predTheta, nav3daction);
        if (cost >= INFINITECOST) {
            continue;
        }

        EnvXXXLATHashEntry_t* OutHashEntry;
        if ((OutHashEntry = (this->*GetHashEntry)(predX, predY, predTheta)) == NULL) {
            // have to create a new entry
            OutHashEntry = (this->*CreateNewHashEntry)(predX, predY, predTheta);
        }

        PredIDV->push_back(OutHashEntry->stateID);
        CostV->push_back(cost);
    }

#if TIME_DEBUG
    time_getsuccs += clock()-currenttime;
#endif
}

void EnvironmentXXXLAT::SetAllActionsandAllOutcomes(CMDPSTATE* state)
{
    int cost;

#if DEBUG
    if (state->StateID >= (int)StateID2CoordTable.size()) {
        throw SBPL_Exception("ERROR in Env... function: stateID illegal");
    }

    if ((int)state->Actions.size() != 0) {
        throw SBPL_Exception("ERROR in Env_setAllActionsandAllOutcomes: actions already exist for the state");
    }
#endif

    // goal state should be absorbing
    if (state->StateID == EnvXXXLAT.goalstateid) {
        return;
    }

    // get X, Y for the state
    EnvXXXLATHashEntry_t* HashEntry = StateID2CoordTable[state->StateID];

    // iterate through actions
    for (int aind = 0; aind < EnvXXXCfg.actionwidth; aind++) {
        EnvXXXLATAction_t* nav3daction = &EnvXXXCfg.ActionsV[(unsigned int)HashEntry->Theta][aind];
        int newX = HashEntry->X + nav3daction->dX;
        int newY = HashEntry->Y + nav3daction->dY;
        int newTheta = normalizeDiscAngle(nav3daction->endtheta);

        // skip the invalid cells
        if (!IsValidCell(newX, newY)) {
            continue;
        }

        // get cost
        cost = GetActionCost(HashEntry->X, HashEntry->Y, HashEntry->Theta, nav3daction);
        if (cost >= INFINITECOST) {
            continue;
        }

        // add the action
        CMDPACTION* action = state->AddAction(aind);

#if TIME_DEBUG
        clock_t currenttime = clock();
#endif

        EnvXXXLATHashEntry_t* OutHashEntry;
        if ((OutHashEntry = (this->*GetHashEntry)(newX, newY, newTheta)) == NULL) {
            // have to create a new entry
            OutHashEntry = (this->*CreateNewHashEntry)(newX, newY, newTheta);
        }
        action->AddOutcome(OutHashEntry->stateID, cost, 1.0);

#if TIME_DEBUG
        time3_addallout += clock()-currenttime;
#endif
    }
}

void EnvironmentXXXLAT::GetPredsofChangedEdges(
    std::vector<nav2dcell_t> const* changedcellsV,
    std::vector<int>* preds_of_changededgesIDV)
{
    nav2dcell_t cell;
    sbpl_xy_theta_cell_t affectedcell;
    EnvXXXLATHashEntry_t* affectedHashEntry;

    // increment iteration for processing savings
    iteration++;

    for (int i = 0; i < (int)changedcellsV->size(); i++) {
        cell = changedcellsV->at(i);

        // now iterate over all states that could potentially be affected
        for (int sind = 0; sind < (int)affectedpredstatesV.size(); sind++) {
            affectedcell = affectedpredstatesV.at(sind);

            // translate to correct for the offset
            affectedcell.x = affectedcell.x + cell.x;
            affectedcell.y = affectedcell.y + cell.y;

            // insert only if it was actually generated
            affectedHashEntry = (this->*GetHashEntry)(affectedcell.x, affectedcell.y, affectedcell.theta);
            if (affectedHashEntry != NULL && affectedHashEntry->iteration < iteration) {
                preds_of_changededgesIDV->push_back(affectedHashEntry->stateID);
                affectedHashEntry->iteration = iteration; // mark as already inserted
            }
        }
    }
}

void EnvironmentXXXLAT::GetSuccsofChangedEdges(
    std::vector<nav2dcell_t> const* changedcellsV,
    std::vector<int>* succs_of_changededgesIDV)
{
    nav2dcell_t cell;
    sbpl_xy_theta_cell_t affectedcell;
    EnvXXXLATHashEntry_t* affectedHashEntry;

    throw SBPL_Exception("ERROR: getsuccs is not supported currently");

    // increment iteration for processing savings
    iteration++;

    // TODO - check
    for (int i = 0; i < (int)changedcellsV->size(); i++) {
        cell = changedcellsV->at(i);

        // now iterate over all states that could potentially be affected
        for (int sind = 0; sind < (int)affectedsuccstatesV.size(); sind++) {
            affectedcell = affectedsuccstatesV.at(sind);

            // translate to correct for the offset
            affectedcell.x = affectedcell.x + cell.x;
            affectedcell.y = affectedcell.y + cell.y;

            // insert only if it was actually generated
            affectedHashEntry = (this->*GetHashEntry)(affectedcell.x, affectedcell.y, affectedcell.theta);
            if (affectedHashEntry != NULL && affectedHashEntry->iteration < iteration) {
                succs_of_changededgesIDV->push_back(affectedHashEntry->stateID);
                // mark as already inserted
                affectedHashEntry->iteration = iteration;
            }
        }
    }
}

void EnvironmentXXXLAT::InitializeEnvironment()
{
    EnvXXXLATHashEntry_t* HashEntry;

    int maxsize = EnvXXXCfg.EnvWidth_c * EnvXXXCfg.EnvHeight_c * EnvXXXCfg.NumThetaDirs;

    if (maxsize <= SBPL_XXXLAT_MAXSTATESFORLOOKUP) {
        SBPL_PRINTF("environment stores states in lookup table\n");

        Coord2StateIDHashTable_lookup = new EnvXXXLATHashEntry_t*[maxsize];
        for (int i = 0; i < maxsize; i++) {
            Coord2StateIDHashTable_lookup[i] = NULL;
        }
        GetHashEntry = &EnvironmentXXXLAT::GetHashEntry_lookup;
        CreateNewHashEntry = &EnvironmentXXXLAT::CreateNewHashEntry_lookup;

        // not using hash table
        HashTableSize = 0;
        Coord2StateIDHashTable = NULL;
    }
    else {
        SBPL_PRINTF("environment stores states in hashtable\n");

        // initialize the map from Coord to StateID
        HashTableSize = 4 * 1024 * 1024; // should be power of two
        Coord2StateIDHashTable = new std::vector<EnvXXXLATHashEntry_t*>[HashTableSize];
        GetHashEntry = &EnvironmentXXXLAT::GetHashEntry_hash;
        CreateNewHashEntry = &EnvironmentXXXLAT::CreateNewHashEntry_hash;

        // not using hash
        Coord2StateIDHashTable_lookup = NULL;
    }

    // initialize the map from StateID to Coord
    StateID2CoordTable.clear();

    // create start state
    if (NULL == (HashEntry = (this->*GetHashEntry)(
            EnvXXXCfg.StartX_c,
            EnvXXXCfg.StartY_c,
            EnvXXXCfg.StartTheta)))
    {
        // have to create a new entry
        HashEntry = (this->*CreateNewHashEntry)(
                EnvXXXCfg.StartX_c,
                EnvXXXCfg.StartY_c,
                EnvXXXCfg.StartTheta);
    }
    EnvXXXLAT.startstateid = HashEntry->stateID;

    // create goal state
    if ((HashEntry = (this->*GetHashEntry)(
            EnvXXXCfg.EndX_c,
            EnvXXXCfg.EndY_c,
            EnvXXXCfg.EndTheta)) == NULL)
    {
        // have to create a new entry
        HashEntry = (this->*CreateNewHashEntry)(
                EnvXXXCfg.EndX_c,
                EnvXXXCfg.EndY_c,
                EnvXXXCfg.EndTheta);
    }
    EnvXXXLAT.goalstateid = HashEntry->stateID;

    // initialized
    EnvXXXLAT.bInitialized = true;
}

// examples of hash functions: map state coordinates onto a hash value
// #define GETHASHBIN(X, Y) (Y*WIDTH_Y+X)
// here we have state coord: <X1, X2, X3, X4>
unsigned int EnvironmentXXXLAT::GETHASHBIN(
    unsigned int X1,
    unsigned int X2,
    unsigned int Theta)
{
    return inthash(inthash(X1) + (inthash(X2) << 1) + (inthash(Theta) << 2)) & (HashTableSize - 1);
}

void EnvironmentXXXLAT::PrintHashTableHist(FILE* fOut)
{
    int s0 = 0, s1 = 0, s50 = 0, s100 = 0, s200 = 0, s300 = 0, slarge = 0;

    for (int j = 0; j < HashTableSize; j++) {
        if ((int)Coord2StateIDHashTable[j].size() == 0)
            s0++;
        else if ((int)Coord2StateIDHashTable[j].size() < 5)
            s1++;
        else if ((int)Coord2StateIDHashTable[j].size() < 25)
            s50++;
        else if ((int)Coord2StateIDHashTable[j].size() < 50)
            s100++;
        else if ((int)Coord2StateIDHashTable[j].size() < 100)
            s200++;
        else if ((int)Coord2StateIDHashTable[j].size() < 400)
            s300++;
        else
            slarge++;
    }
    SBPL_FPRINTF(fOut, "hash table histogram: 0:%d, <5:%d, <25:%d, <50:%d, <100:%d, <400:%d, >400:%d\n", s0, s1, s50,
                 s100, s200, s300, slarge);
}

int EnvironmentXXXLAT::GetFromToHeuristic(int FromStateID, int ToStateID)
{
#if USE_HEUR==0
    return 0;
#endif

#if DEBUG
    if (FromStateID >= (int)StateID2CoordTable.size() ||
        ToStateID >= (int)StateID2CoordTable.size())
    {
        SBPL_ERROR("ERROR in EnvXXXLAT... function: stateID illegal\n");
        throw SBPL_Exception();
    }
#endif

    // get X, Y for the state
    EnvXXXLATHashEntry_t* FromHashEntry = StateID2CoordTable[FromStateID];
    EnvXXXLATHashEntry_t* ToHashEntry = StateID2CoordTable[ToStateID];

    // TODO - check if one of the gridsearches already computed and then use it.

    return (int)(XXXLAT_COSTMULT_MTOMM *
            EuclideanDistance_m(FromHashEntry->X, FromHashEntry->Y, ToHashEntry->X, ToHashEntry->Y) /
            EnvXXXCfg.nominalvel_mpersecs);

}

int EnvironmentXXXLAT::GetGoalHeuristic(int stateID)
{
#if USE_HEUR==0
    return 0;
#endif

#if DEBUG
    if (stateID >= (int)StateID2CoordTable.size()) {
        throw SBPL_Exception("ERROR in EnvXXXLAT... function: stateID illegal");
    }
#endif

    EnvXXXLATHashEntry_t* HashEntry = StateID2CoordTable[stateID];
    // computes distances from start state that is grid2D, so it is EndX_c EndY_c
    int h2D = grid2Dsearchfromgoal->getlowerboundoncostfromstart_inmm(HashEntry->X, HashEntry->Y);
    int hEuclid = (int)(XXXLAT_COSTMULT_MTOMM *
            EuclideanDistance_m(HashEntry->X, HashEntry->Y, EnvXXXCfg.EndX_c, EnvXXXCfg.EndY_c));

    // Add clearance component based on Grid2D cost
    double clearancePenalty = 0;
    const int CLEARANCE_WINDOW = 50;  // Increased window to look further for obstacles
    double minCostRatio = 1.0;  // Track the best clearance we find
    // First find minimum cost ratio (best clearance) in the window
    for(int dx = -CLEARANCE_WINDOW; dx <= CLEARANCE_WINDOW; dx++) {
        for(int dy = -CLEARANCE_WINDOW; dy <= CLEARANCE_WINDOW; dy++) {
            int checkX = HashEntry->X + dx;
            int checkY = HashEntry->Y + dy;
            if(IsValidCell(checkX, checkY)) {
                int cellCost = EnvXXXCfg.Grid2D[checkX][checkY];
                double costRatio = (double)cellCost / EnvXXXCfg.cost_inscribed_thresh;
                double distance = sqrt(dx*dx + dy*dy);
                // For each direction, accumulate inverse clearance
                clearancePenalty += costRatio * (CLEARANCE_WINDOW - distance) / CLEARANCE_WINDOW;
            }
        }
    }
    // Scale the clearance penalty - even small costs will contribute
    const double CLEARANCE_WEIGHT = 1.0;  
    int totalHeuristic = (int)(((double)__max(h2D, hEuclid)) / EnvXXXCfg.nominalvel_mpersecs + 
                              CLEARANCE_WEIGHT * clearancePenalty);
    return totalHeuristic;

    // define this function if it is used in the planner (heuristic backward search would use it)
    // return (int)(((double)__max(h2D, hEuclid)) / EnvXXXCfg.nominalvel_mpersecs);
}

int EnvironmentXXXLAT::GetStartHeuristic(int stateID)
{
#if USE_HEUR==0
    return 0;
#endif

#if DEBUG
    if (stateID >= (int)StateID2CoordTable.size()) {
        throw SBPL_Exception("ERROR in EnvXXXLAT... function: stateID illegal");
    }
#endif

    EnvXXXLATHashEntry_t* HashEntry = StateID2CoordTable[stateID];
    int h2D = grid2Dsearchfromstart->getlowerboundoncostfromstart_inmm(HashEntry->X, HashEntry->Y);
    int hEuclid = (int)(XXXLAT_COSTMULT_MTOMM *
            EuclideanDistance_m(EnvXXXCfg.StartX_c, EnvXXXCfg.StartY_c, HashEntry->X, HashEntry->Y));

    // define this function if it is used in the planner (heuristic backward search would use it)
    return (int)(((double)__max(h2D, hEuclid)) / EnvXXXCfg.nominalvel_mpersecs);
}

int EnvironmentXXXLAT::SizeofCreatedEnv()
{
    return (int)StateID2CoordTable.size();
}

const EnvXXXLATHashEntry_t*
EnvironmentXXXLAT::GetStateEntry(int state_id) const
{
    if (state_id >= 0 && state_id < (int)StateID2CoordTable.size()) {
        return StateID2CoordTable[state_id];
    }
    else {
        return NULL;
    }
}

//------------------------------------------------------------------------------


void EnvironmentXXXLAT::GetLazySuccs(
    int SourceStateID,
    std::vector<int>* SuccIDV,
    std::vector<int>* CostV,
    std::vector<bool>* isTrueCost,
    std::vector<EnvXXXLATAction_t*>* actionV)
{
    int aind;

#if TIME_DEBUG
    clock_t currenttime = clock();
#endif

    // clear the successor array
    SuccIDV->clear();
    CostV->clear();
    SuccIDV->reserve(EnvXXXCfg.actionwidth);
    CostV->reserve(EnvXXXCfg.actionwidth);
    isTrueCost->reserve(EnvXXXCfg.actionwidth);
    if (actionV != NULL) {
        actionV->clear();
        actionV->reserve(EnvXXXCfg.actionwidth);
    }

    // goal state should be absorbing
    if (SourceStateID == EnvXXXLAT.goalstateid) {
        return;
    }

    // get X, Y for the state
    EnvXXXLATHashEntry_t* HashEntry = StateID2CoordTable[SourceStateID];

    // iterate through actions
    for (aind = 0; aind < EnvXXXCfg.actionwidth; aind++) {
        EnvXXXLATAction_t* nav3daction = &EnvXXXCfg.ActionsV[(unsigned int)HashEntry->Theta][aind];
        int newX = HashEntry->X + nav3daction->dX;
        int newY = HashEntry->Y + nav3daction->dY;
        int newTheta = normalizeDiscAngle(nav3daction->endtheta);

        // skip the invalid cells
        if (!IsValidCell(newX, newY)) {
            continue;
        }

        // if we are supposed to return the action, then don't do lazy
        if (!actionV) {
            EnvXXXLATHashEntry_t* OutHashEntry;
            if ((OutHashEntry = (this->*GetHashEntry)(newX, newY, newTheta)) == NULL) {
                OutHashEntry = (this->*CreateNewHashEntry)(newX, newY, newTheta);
            }
            SuccIDV->push_back(OutHashEntry->stateID);
            CostV->push_back(nav3daction->cost);
            isTrueCost->push_back(false);
            continue;
        }

        // get cost
        int cost = GetActionCost(HashEntry->X, HashEntry->Y, HashEntry->Theta, nav3daction);
        if (cost >= INFINITECOST) {
            continue;
        }

        EnvXXXLATHashEntry_t* OutHashEntry;
        if ((OutHashEntry = (this->*GetHashEntry)(newX, newY, newTheta)) == NULL) {
            // have to create a new entry
            OutHashEntry = (this->*CreateNewHashEntry)(newX, newY, newTheta);
        }

        SuccIDV->push_back(OutHashEntry->stateID);
        CostV->push_back(cost);
        isTrueCost->push_back(true);
        if (actionV != NULL) {
            actionV->push_back(nav3daction);
        }
    }

#if TIME_DEBUG
    time_getsuccs += clock()-currenttime;
#endif
}

int EnvironmentXXXLAT::GetTrueCost(int parentID, int childID)
{
    EnvXXXLATHashEntry_t* fromHash = StateID2CoordTable[parentID];
    EnvXXXLATHashEntry_t* toHash = StateID2CoordTable[childID];

    for (int i = 0; i < EnvXXXCfg.actionwidth; i++) {
        EnvXXXLATAction_t* nav3daction = &EnvXXXCfg.ActionsV[(unsigned int)fromHash->Theta][i];
        int newX = fromHash->X + nav3daction->dX;
        int newY = fromHash->Y + nav3daction->dY;
        int newTheta = normalizeDiscAngle(nav3daction->endtheta);

        // skip the invalid cells
        if (!IsValidCell(newX, newY)) {
            continue;
        }

        EnvXXXLATHashEntry_t* hash;
        if ((hash = (this->*GetHashEntry)(newX, newY, newTheta)) == NULL) {
            continue;
        }
        if (hash->stateID != toHash->stateID) {
            continue;
        }

        // get cost
        int cost = GetActionCost(fromHash->X, fromHash->Y, fromHash->Theta, nav3daction);

        if (cost >= INFINITECOST) {
            return -1;
        }
        return cost;
    }
    throw SBPL_Exception("This should never happen! we didn't find the state we need to get the true cost for!");
    return -1;
}

void EnvironmentXXXLAT::GetSuccsWithUniqueIds(
    int SourceStateID,
    std::vector<int>* SuccIDV,
    std::vector<int>* CostV,
    std::vector<EnvXXXLATAction_t*>* actionV)
{
    GetSuccs(SourceStateID, SuccIDV, CostV, actionV);
}

void EnvironmentXXXLAT::GetLazySuccsWithUniqueIds(
    int SourceStateID,
    std::vector<int>* SuccIDV,
    std::vector<int>* CostV,
    std::vector<bool>* isTrueCost,
    std::vector<EnvXXXLATAction_t*>* actionV)
{
    GetLazySuccs(SourceStateID, SuccIDV, CostV, isTrueCost, actionV);
}

bool EnvironmentXXXLAT::isGoal(int id)
{
    return EnvXXXLAT.goalstateid == id;
}

void EnvironmentXXXLAT::GetLazyPreds(
    int TargetStateID,
    std::vector<int>* PredIDV,
    std::vector<int>* CostV,
    std::vector<bool>* isTrueCost)
{
    int aind;

#if TIME_DEBUG
    clock_t currenttime = clock();
#endif

    // get X, Y for the state
    EnvXXXLATHashEntry_t* HashEntry = StateID2CoordTable[TargetStateID];

    // clear the successor array
    PredIDV->clear();
    CostV->clear();
    PredIDV->reserve(EnvXXXCfg.PredActionsV[(unsigned int)HashEntry->Theta].size());
    CostV->reserve(EnvXXXCfg.PredActionsV[(unsigned int)HashEntry->Theta].size());

    // iterate through actions
    std::vector<EnvXXXLATAction_t*>* actionsV = &EnvXXXCfg.PredActionsV[(unsigned int)HashEntry->Theta];
    for (aind = 0; aind < (int)EnvXXXCfg.PredActionsV[(unsigned int)HashEntry->Theta].size(); aind++)
    {
        EnvXXXLATAction_t* nav3daction = actionsV->at(aind);

        int predX = HashEntry->X - nav3daction->dX;
        int predY = HashEntry->Y - nav3daction->dY;
        int predTheta = nav3daction->starttheta;

        //skip the invalid cells
        if (!IsValidCell(predX, predY)) {
            continue;
        }

        EnvXXXLATHashEntry_t* OutHashEntry;
        if ((OutHashEntry = (this->*GetHashEntry)(predX, predY, predTheta)) == NULL) {
            OutHashEntry = (this->*CreateNewHashEntry)(predX, predY, predTheta);
        }

        PredIDV->push_back(OutHashEntry->stateID);
        CostV->push_back(nav3daction->cost);
        isTrueCost->push_back(false);
    }

#if TIME_DEBUG
    time_getsuccs += clock()-currenttime;
#endif
}

void EnvironmentXXXLAT::GetPredsWithUniqueIds(
    int TargetStateID,
    std::vector<int>* PredIDV,
    std::vector<int>* CostV)
{
    GetPreds(TargetStateID, PredIDV, CostV);
}

void EnvironmentXXXLAT::GetLazyPredsWithUniqueIds(
    int TargetStateID,
    std::vector<int>* PredIDV,
    std::vector<int>* CostV,
    std::vector<bool>* isTrueCost)
{
    GetLazyPreds(TargetStateID, PredIDV, CostV, isTrueCost);
}
