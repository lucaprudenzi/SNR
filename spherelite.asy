settings.outformat = "pdf";
settings.prc=false;
settings.render=16;
size(300);

import three;
// In this program there is the following convetion:
// X axis : blue axis
// Y axis : green axis
// Z axis : red axis

///////////////////////////////////////////////////////////////////////
// PARAMETERS IN GEOCENTRIC FRAME
file geoparamsfile = input("geoparams.txt");
real [] geoparams = geoparamsfile;

real thetaS = geoparams[4]*pi/180;
real phiS = geoparams[5]*pi/180;
real psi = geoparams[6]*pi/180;
real iota = geoparams[7]*pi/180;
// 1 for H1 detector, 2 for L1 detector, 3 for V1 detector
real detector = geoparams[8];
real theta;
real phi;
real south;

// Hanford coordinates
if (detector == 1){
	real latitude = 46+27/60.+19/3600;
	real longitude = 2*pi-(119+44/60+28/3600);
	// x-arm orientation with respect to South position
	real rotation = 180+36;
	
	theta = (90-latitude)*pi/180;
	phi = longitude*pi/180;
	south = rotation*pi/180;
}
// Livingstone coordinates
if (detector == 2){
	real latitude = 30+33/60+46/3600;
	real longitude = 2*pi-(90+46/60+27/3600);
	// x-arm orientation with respect to South position
	real rotation = 270+18;

	theta = (90-latitude)*pi/180;
	phi = longitude*pi/180;
	south = rotation*pi/180;
}

// Virgo coordinates
if (detector == 3){
	real latitude = 43+37/60+53/3600;
	real longitude = 10+30/60+16/3600;
	// x-arm orientation with respect to South position
	real rotation = 180-19;

	theta = (90-latitude)*pi/180;
	phi = longitude*pi/180;
	south = rotation*pi/180;
}

// Test  variables
// Detector coordinates
//theta=30*pi/180;//90-Latitude
//phi=30*pi/180;//Longitude
//south=10*pi/180;//rotation angle of x-axis from south

// Source coordinates 
//real thetaS=45*pi/180;//90-Latitude
//real phiS=50*pi/180;//Longitude

// binary plane in the sky 
//real iota = 30*pi/180; // theta for binary (direction of z-ax respect to z-geocentric)
//real psi = 80*pi/180; // phi for binary (between x-ax binary and projected x-ax of geocentric frame)

///////////////////////////////////////////////////////////////////////

// geocentric axes
draw((0,0,0)--(1.5,0,0), blue, arrow=Arrow3()); //x-axis
draw((0,0,0)--(0,1.5,0), green, arrow=Arrow3()); //y-axis
draw((0,0,0)--(0,0,1.5), red, arrow=Arrow3()); //z-axis

// Z axis direction  after second rotation (axis for third rotation)
triple Zdot = rotate(theta*180/pi,(-sin(phi),cos(phi),0))*Z;
// Y axis direction after first rotation (axis for second rotation)
triple Ydot = rotate(phi*180/pi,Z)*Y;

// Detector plane
path3 detplane=(0.5,0,0)--(0,0.5,0)--(-0.5,0,0)--(0,-0.5,0)--cycle;
detplane = rotate(phi*180/pi,Z)*detplane;
detplane = rotate(theta*180/pi,Ydot)*detplane;
detplane = rotate(south*180/pi,Zdot)*detplane;
draw(detplane, 0.1cyan);

// Detector axes
path3 Xax = (0,0,0)--(0.5,0,0);
Xax = rotate(phi*180/pi,Z)*Xax;
Xax = rotate(theta*180/pi,Ydot)*Xax;
Xax = rotate(south*180/pi,Zdot)*Xax;

path3 Yax = (0,0,0)--(0,0.5,0);
Yax = rotate(phi*180/pi,Z)*Yax;
Yax = rotate(south*180/pi,Zdot)*Yax;

path3 Zax = (0,0,0)--(0,0,0.5);
Zax = rotate(theta*180/pi,Ydot)*Zax;

draw(Xax, arrow=Arrow3(), blue);
draw(Yax, arrow=Arrow3(), green);
draw(Zax, arrow=Arrow3(), red);

// Sphere
draw(scale3(1.1)*unitsphere, opacity(0.1));

// Source position
real radius = 1;
triple Source=radius*(sin(thetaS)*cos(phiS),sin(thetaS)*sin(phiS),cos(thetaS));
triple projectedSource=radius*(sin(thetaS)*cos(phiS),sin(thetaS)*sin(phiS),0);
dot(Source);
label("S",Source,N);
dot(projectedSource);

// source projection on geocentric system
path3 projectedSourceXY = O--projectedSource;
path3 projectedSourceZ = Source--projectedSource;
draw(projectedSourceZ);
draw(projectedSourceXY);
draw(O--Source);

path3 anglearc(real radius, triple A, triple B, triple C) {
   triple center = B;
   triple start = B + radius * unit(A-B);
   return arc(center, start, C);
}

// Source Angles on geocentric system
path3 arcphiS = anglearc(0.5, X, O, projectedSource);
path3 arcthetaS = anglearc(0.5, Z, O, Source);

draw(arcphiS, L=Label("$\small\phi_s$", align=N));
draw(arcthetaS, L=Label("$\small\theta_s$", align=N));

// Two methods to find out theta and phi in detector frame
/*
// METHOD 1 (INTERSECTION METHOD) for angles in detector frame

// Surface identical to detplane, bigger to be sure of intersection
// ghostsurface is used for intersection with intersectionpoints(path3,surface)
path3 ghostplane=(2.5,0,0)--(0,2.5,0)--(-2.5,0,0)--(0,-2.5,0)--cycle;
ghostplane = rotate(phi*180/pi,Z)*ghostplane;
ghostplane = rotate(theta*180/pi,Ydot)*ghostplane;
ghostplane = rotate(south*180/pi,Zdot)*ghostplane;
surface ghostsurface=surface(ghostplane);
//draw(ghostsurface);

// Normal to the new plane
real newZaxlenght = 5; // a sufficient long line
path3 newZax = (0,0,+newZaxlenght)--(0,0,-newZaxlenght); // -newZaxlenght because shift moves head
newZax = rotate(theta*180/pi,Ydot)*newZax;
newZax = shift(Source)*newZax;
//draw(newZax);

// intersection point between normal to detector plane and detector xy plane
triple [] ip = intersectionpoints(newZax, ghostsurface);
path3 proiezione1 = Source--ip[0];
path3 proiezione2 = ip[0]--O;
draw(proiezione1);
draw(proiezione2);
dot(ip[0]);
//Position of X axis after three rotation
triple Xdot = rotate(phi*180/pi,Z)*X;
Xdot = rotate(theta*180/pi,Ydot)*Xdot;
Xdot = rotate(south*180/pi,Zdot)*Xdot;

// Angles on detector plane
path3 arcphiSD = anglearc(0.2, Xdot, O, ip[0]);
path3 arcthetaSD = anglearc(0.2, Zdot, O, Source);
draw(arcphiSD);
draw(arcthetaSD);
label("$\tiny\phi_{d}$",arcphiSD,N);
label("$\tiny\theta_{d}$",arcthetaSD,N);

triple detSource = rotate(-south*180/pi,Zdot)*Source;
detSource = rotate(-theta*180/pi,Ydot)*detSource;
detSource = rotate(-phi*180/pi,Z)*detSource;
real longSource = longitude(detSource);
real latSource = latitude(detSource);
// Print lat and long of the Source in detector frame
write("latitude in detector frame: ", latSource);
write("longitude in detector frame: ", longSource);

*/

// METHOD 2 (DIRECT COMPUTATION) for angles in detector frame
// Print angles of detector frame 
triple detSource = rotate(-south*180/pi,Zdot)*Source;
detSource = rotate(-theta*180/pi,Ydot)*detSource;
detSource = rotate(-phi*180/pi,Z)*detSource;
real longSource = longitude(detSource);
real latSource = latitude(detSource);
// Print lat and long of the Source in detector frame
write(90-latSource); //theta_d
write(longSource); //phi_d

real thetaSD = (90-latSource)*pi/180;
real phiSD = longSource*pi/180; 
triple detSource=radius*(sin(thetaSD)*cos(phiSD),sin(thetaSD)*sin(phiSD),cos(thetaSD));
triple detprojectedSource=radius*(sin(thetaSD)*cos(phiSD),sin(thetaSD)*sin(phiSD),0);

path3 arcphiSD = anglearc(0.2, X, O, detprojectedSource);
path3 arcthetaSD = anglearc(0.2, Z, O, detSource);

arcphiSD = rotate(phi*180/pi,Z)*arcphiSD;
arcphiSD = rotate(theta*180/pi,Ydot)*arcphiSD;
arcphiSD = rotate(south*180/pi,Zdot)*arcphiSD;

arcthetaSD = rotate(phi*180/pi,Z)*arcthetaSD;
arcthetaSD = rotate(theta*180/pi,Ydot)*arcthetaSD;
arcthetaSD = rotate(south*180/pi,Zdot)*arcthetaSD;

draw(arcphiSD, L=Label("$\small\phi_d$", align=N));
draw(arcthetaSD, L=Label("$\small\theta_d$", align=N));

detprojectedSource = rotate(phi*180/pi,Z)*detprojectedSource;
detprojectedSource = rotate(theta*180/pi,Ydot)*detprojectedSource;
detprojectedSource = rotate(south*180/pi,Zdot)*detprojectedSource;

dot(detprojectedSource);
path3 detprojectedSourceXY = O--detprojectedSource;
path3 detprojectedSourceZ = Source--detprojectedSource;
draw(detprojectedSourceZ);
draw(detprojectedSourceXY);

// BINARY PLANE SECTION
// We have to find out what is phi angle for binary system knowing psi

// PART1: Source plane and geocentric plane
// Binary plane
path3 binaryplane=(0.5,0,0)--(0,0.5,0)--(-0.5,0,0)--(0,-0.5,0)--cycle;
// Binary axes
path3 binXax = (0,0,0)--(0.5,0,0);
path3 binYax = (0,0,0)--(0,0.5,0);
path3 binZax = (0,0,0)--(0,0,0.5);

//point for direct rotation of binary plane
triple projectedSourceRot = rotate(90,(0,0,1))*projectedSource;

// costruction of a surface equal to binary plane, used for intersection
path3 binAuxiliarplane=(2.5,0,0)--(0,2.5,0)--(-2.5,0,0)--(0,-2.5,0)--cycle;
binAuxiliarplane = rotate((thetaS+iota)*180/pi,projectedSourceRot)*binAuxiliarplane;
surface binAuxiliarsurface=surface(binAuxiliarplane);

// Construction of normal to binary plane passing by X axis of geocentric plane
real binprojectionlenght = 2; // a sufficient long line
path3 binprojection = (0,0,binprojectionlenght)--(0,0,-binprojectionlenght);
binprojection = rotate((thetaS+iota)*180/pi,(-sin(phiS),cos(phiS),0))*binprojection;
binprojection = shift(0.5,0,0)*binprojection;
//draw(binprojection);

// Intersection point between binary plane and the normal passing through x axis of geocentric plane
triple [] binprojecteddot = intersectionpoints(binprojection, binAuxiliarsurface);
//dot(binprojecteddot[0]);
binprojecteddot[0] = rotate(-(thetaS+iota)*180/pi,projectedSourceRot)*binprojecteddot[0];
//dot(binprojecteddot[0]);
// psi geocentric
real binlong = longitude(binprojecteddot[0])*pi/180;

triple binposition = rotate(psi*180/pi,(0,0,1))*binprojecteddot[0];

// binphi is the angle of rotation of binary xy plane
real binphi = longitude(binposition);
binphi = binphi*pi/180;

// Now I orient binary plane in the right way
binaryplane = rotate(binphi*180/pi,Z)*binaryplane;
binaryplane = rotate((thetaS+iota)*180/pi,(-sin(phiS),cos(phiS),0))*binaryplane;
binXax = rotate(binphi*180/pi,Z)*binXax;
binXax = rotate((thetaS+iota)*180/pi,(-sin(phiS),cos(phiS),0))*binXax;
binYax = rotate(binphi*180/pi,Z)*binYax;
binYax = rotate((thetaS+iota)*180/pi,(-sin(phiS),cos(phiS),0))*binYax;
binZax = rotate(binphi*180/pi,Z)*binZax;
binZax = rotate((thetaS+iota)*180/pi,(-sin(phiS),cos(phiS),0))*binZax;

//psi arc (geocentric frame)
path3 arcpsi = anglearc(0.3, binprojecteddot[0],O, binposition);
arcpsi = rotate((thetaS+iota)*180/pi,(-sin(phiS),cos(phiS),0))*arcpsi;
triple binprojectionlinedot = rotate((thetaS+iota)*180/pi,(-sin(phiS),cos(phiS),0))*binprojecteddot[0];
path3 binnormal = O--binprojectionlinedot;

// iota arc
path3 arciota = anglearc(0.2,Source, O, rotate((thetaS+iota)*180/pi,(-sin(phiS),cos(phiS),0))*Z);

// PART2: Source plane and detector plane

// X-axis of detector after detector rotation 
triple detXdot = rotate(phi*180/pi,Z)*(0.5,0,0);
detXdot = rotate(theta*180/pi,Ydot)*detXdot;
detXdot = rotate(south*180/pi,Zdot)*detXdot;
// X-axis of binary after binary rotation
triple binXdot = rotate(binphi*180/pi,Z)*(0.5,0,0);
binXdot = rotate((thetaS+iota)*180/pi,(-sin(phiS),cos(phiS),0))*binXdot;

// Shift the normal to binary plane to Xaxis of detector system
path3 binprojection = (0,0,binprojectionlenght)--(0,0,-binprojectionlenght);
binprojection = rotate((thetaS+iota)*180/pi,(-sin(phiS),cos(phiS),0))*binprojection;
binprojection = shift(detXdot)*binprojection;
//draw(binprojection);

// Intersection point between binary plane and the normal passing through x axis of detector
triple [] bindetprojecteddot = intersectionpoints(binprojection, binAuxiliarsurface);
//dot(bindetprojecteddot[0]);
// Rotation of intersection point in geocentric plane
bindetprojecteddot[0] = rotate(-(thetaS+iota)*180/pi,projectedSourceRot)*bindetprojecteddot[0];
bindetprojecteddot[0] = rotate(-binphi*180/pi,Z)*bindetprojecteddot[0];
//dot(bindetprojecteddot[0]);
// psi_det
real bindetlong = longitude(bindetprojecteddot[0]);
write(360-bindetlong); //psi

//psidet arc
path3 psidetarc = anglearc(0.2,bindetprojecteddot[0],O,X);
path3 bindetnormal = (O--bindetprojecteddot[0]);
psidetarc = rotate(binphi*180/pi,Z)*psidetarc;
psidetarc = rotate((thetaS+iota)*180/pi,projectedSourceRot)*psidetarc;
bindetnormal = rotate(binphi*180/pi,Z)*bindetnormal;
bindetnormal = rotate((thetaS+iota)*180/pi,(-sin(phiS),cos(phiS),0))*bindetnormal;

// Shifting all objects to Source position
binaryplane = shift(Source)*binaryplane;
binXax = shift(Source)*binXax;
binYax = shift(Source)*binYax;
binZax = shift(Source)*binZax;
psidetarc = shift(Source)*psidetarc;
arciota = shift(Source)*arciota;
arcpsi = shift(Source)*arcpsi;
bindetnormal = shift(Source)*bindetnormal;
binnormal = shift(Source)*binnormal;

// normal to binary plane (used for draw iota angle)
path3 binspherenormal = O--(0,0,0.5);
binspherenormal = rotate(thetaS*180/pi,(-sin(phiS),cos(phiS),0))*binspherenormal;
binspherenormal = shift(Source)*binspherenormal;
draw(binspherenormal,arrow=Arrow3(),red);

draw(psidetarc);
label("$\psi_d$",psidetarc,N);
draw(arcpsi);
label("$\psi_s$",arcpsi,N);
draw(arciota);
label("$\iota$",arciota,N);
draw(binnormal);
draw(bindetnormal);
draw(binaryplane);
draw(binXax, arrow=Arrow3(), blue);
draw(binYax, arrow=Arrow3(), green);
draw(binZax, arrow=Arrow3(), red);


