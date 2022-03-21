
mm = 1e-3;
R = 1*mm;
R1 = 0.4 * mm;
R2 = 0.35 * mm;
R3 = 0.25 * mm;

L = 25*mm;
L1 = 0.125*mm;
L2 = 2*mm;
L3 = 2.25*mm;
rr = 0.5*mm;

lc = 0.1*mm;

Point(1) = {L1,0,0,lc/3};
Point(2) = {L,0,0,2*lc};
Point(3) = {L,R3,0,2*lc};
Point(4) = {L3,R3,0,lc};
Point(5) = {L2,R2,0,lc/2};
Point(6) = {L1,R2,0,lc/4};
Point(7) = {L1,R1,0,lc/4};
Point(8) = {0, R1, 0, lc/4};
Point(9) = {0, R,  0, lc};
Point(10) = {L2-rr,R,0,lc};
Point(11) = {L2-rr,R-rr, 0, lc};
Point(12) = {L2, R-rr, 0, lc};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,1};

Line(7) = {6,7};
Line(8) = {7,8};
Line(9) = {8,9};
Line(10) = {9,10};
Circle(11) = {10,11,12};
Line(12) = {12,5};

Curve Loop(21) = {1,2,3,4,5,6};
Plane Surface(22) = {21};
Curve Loop(23) = {5,7,8,9,10,11,12};
Plane Surface(24) = {23};

Physical Curve("hsteel") = {2,3,4};
Physical Curve("temp_steel") = {6};
Physical Curve("temp_glass") = {7,8};
Physical Curve("hglass") = {10,11,12};

Physical Surface("steel") = {22};
Physical Surface("glass") = {24};







