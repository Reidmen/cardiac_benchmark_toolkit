Mesh.CharacteristicLengthMin = 3;
Mesh.CharacteristicLengthMax = 3;
Merge "left_ventricle.stl";
Merge "left_ventricle_lit.stl";
Merge "right_ventricle.stl";
Merge "right_ventricle_lit.stl";
Coherence Mesh;
Surface Loop(1) = {1, 2};
Surface Loop(2) = {3, 4};
Volume(1) = {1};
Volume(2) = {2};
//+
Physical Surface("left_ventricle", 5) = {1};
//+
Physical Surface("left_ventricle_lit", 6) = {2};
//+
Physical Surface("right_ventricle", 7) = {3};
//+
Physical Surface("right_ventricle_lit", 8) = {3};
//+
Physical Volume("left_right_lumen", 9) = {1, 2};
