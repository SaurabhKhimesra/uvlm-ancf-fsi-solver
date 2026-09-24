%% path

addpath ./save;
addpath ./cores;
addpath ./cores/functions;
addpath ./cores/functions/structure;
addpath ./cores/functions/fluid;
addpath ./cores/solver;
addpath ./cores/solver/structure;
addpath ./cores/solver/fluid;

%% ToolBox

addpath ./cores/ToolBoxes/Plate_Mesh;
% addpath ./cores/ToolBoxes/mpg_write/src;                % only for movie_format = 'mpeg'; get_toolboxes does not fetch it
addpath ./cores/ToolBoxes/mmwrite;
% addpath ./cores/ToolBoxes/lightspeed;
addpath ./cores/ToolBoxes/mntimes;
addpath ./cores/ToolBoxes/TriStream;
addpath ./cores/ToolBoxes/quiver5;
addpath ./cores/ToolBoxes/FEM_sparse/FEM_sparse;