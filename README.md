sbpl_towing
===========================

This is a modified version of the original `sbpl` (https://github.com/sbpl/sbpl)

## Build

I recommend creating a new catkin workspace for this repository and the ROS wrapper repo (https://github.com/jordi-lete/sbpl_lattice_planner_towing) 
if you haven't done this already:

```bash
mkdir -p ~/sbpl_towing_ws/src
cd ~/sbpl_towing_ws/
catkin_make
```

Then clone this repo and build

```bash
cd src/
git clone https://github.com/jordi-lete/sbpl_towing.git
```

```bash
cd sbpl_towing/
mkdir build
cd build/
cmake ..
make
cd ..
source devel/setup.bash
```

You can also add the source command to your ./bashrc script to ensure the workspace is sourced for all future shells.

## SBPL documentation:

I. Building, Installing, and Using SBPL

    SBPL is available as a standalone software library. SBPL itself has no
    dependencies other than the C/C++ standard library.

    These build and install instructions are primarily for Linux. For other
    operating systems, CMake can generate the platform-specific build and project
    files necessary for building SBPL.

    Versions of ROS older than Fuerte may contain packages that depend on a ROS
    package version of SBPL. The recommended method to install SBPL is to install
    it as a standard system library. However, if you wish to use the old ROS
    package version of SBPL, you may follow these instructions.

    1. Building and Installing SBPL from source

        1.1 Build SBPL

            SBPL uses git as its version control system. From the directory where
            you want the SBPL source to reside, clone the latest source from
            https://wmg-gitlab.wmgds.wmg.warwick.ac.uk/intelligent-vehicles/sbpl_towing:

            git clone git@wmg-gitlab.wmgds.wmg.warwick.ac.uk:intelligent-vehicles/sbpl_towing.git

            In the source directory, build the SBPL library using standard
            CMake build conventions:

            mkdir build
            cd build
            cmake .. -DCMAKE_INSTALL_PREFIX=$(pwd)/../install
            make
            make install

II. Usage

    Examples for how to use SBPL are in src/test/main.cpp.  Please follow the
    examples carefully. The library contains a number of planning problem
    examples, stored as ascii files.  These files (with extension .cfg) should
    be passed in as arguments into the main function in main.cpp. The files
    can be found in env_examples directory.

    Command-line usage for the test_sbpl program can be viewed by passing '-h'
    as argument to the executable.

    Examples:

        The following can be run from the directory containing test_sbpl,
        which we assume is a build directory in the root of this project.

        $ ./test_sbpl ../env_examples/nav3d/env1.cfg
        Environment: xytheta; Planner: arastar; Search direction: backward
        Initializing ARAPlanner...
        start planning...
        done planning
        size of solution=16
        solution size=0
        Solution is found

        $ ./test_sbpl --env=2d ../env_examples/nav2d/env1.cfg #2d is needed here in order to use 2d config
        Environment: 2d; Planner: arastar; Search direction: backward
        Initializing ARAPlanner...
        start planning...
        done planning
        size of solution=22
        Solution is found

        $ ./test_sbpl --env=robarm --search-dir=forward --planner=rstar ../env_examples/robarm/env1_6d.cfg
        Environment: robarm; Planner: rstar; Search direction: forward
        Initializing RSTARPlanner...
        start planning...
        done planning
        size of solution=44
        Solution is found

   Motion primitives files can be found in sbpl/matlab/mprim directory.

    Finally, few visualization scripts can be found in
    sbpl/matlab/visualization. In particular, plot_3Dpath.m function can be
    used to visualize the path found by xytheta lattice planner. This
    functions takes in .cfg file that specified environment and sol.txt file
    that was generated within main.cpp by xythetalattice planners.

    Note: If you compile the library with the ROS symbol defined, all text
    output will be redirected to ROS logging constructions. Without the ROS
    symbol defined, SBPL will print messages to stdout and test_sbpl will
    generate a solution file, sol.txt, as well as a debugging information
    file, debug.txt

III. Links

    These instructions and more tutorials can be found at www.sbpl.net

    For more information and documentation on SBPL visit:

    http://www.ros.org/wiki/sbpl

    For more information and documentation on using the x,y,theta environment
    available in ROS visit:

    http://www.ros.org/wiki/sbpl_lattice_planner
