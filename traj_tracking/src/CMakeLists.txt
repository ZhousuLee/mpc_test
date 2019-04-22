INCLUDE_DIRECTORIES(/home/lizhousu/桌面/mpc_test/traj_tracking/include  /usr/include/eigen3  /usr/include/)
LINK_DIRECTORIES(/usr/local/lib)
ADD_DEFINITIONS(-std=c++11)
ADD_EXECUTABLE(main main.cpp traj_generate.cpp mpc_control.cpp sim_locate.cpp)
TARGET_LINK_LIBRARIES(main libqpOASES.a)

ADD_EXECUTABLE(sim_locate_test sim_locate_test.cpp traj_generate.cpp mpc_control.cpp sim_locate.cpp)
TARGET_LINK_LIBRARIES(sim_locate_test libqpOASES.a)
