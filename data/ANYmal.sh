cp -r ./ANYmal ../../diffne-build
cd ../../diffne-build
./mainArticulatedLoader ANYmal visualMesh=1
./mainSimplifiedDynamics ANYmal symDir=1 zRange=0,0,0.025 zMargin=3.1
cd ../diffne


