cp -r ./LittleDog ../../diffne-build
cd ../../diffne-build
./mainArticulatedLoader LittleDog visualMesh=1
./mainSimplifiedDynamics LittleDog symDir=1 zRange=0,0,0.02 zMargin=2
cd ../diffne
