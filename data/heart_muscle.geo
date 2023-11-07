Mesh.CharacteristicLengthMin = 2;
Mesh.CharacteristicLengthMax = 2;
Merge "left_endocardium.stl";
Merge "right_endocardium.stl";
Merge "left_epicardium.stl";
Merge "right_epicardium.stl";
Merge "base_and_septum_epicardium.stl";
Coherence Mesh;
Surface Loop(1) = {1, 2, 3, 4, 5};
Volume(1) = {1};
//+
Physical Surface("right_endocardium", 6) = {2};
//+
Physical Surface("left_endocardium", 7) = {1};
//+
Physical Surface("right_epicardium", 8) = {4};
//+
Physical Surface("baser_and_septum_epicardium", 9) = {5};
//+
Physical Surface("left_epicardium", 10) = {3};
//+
Physical Volume("heart_muscle", 11) = {1};
