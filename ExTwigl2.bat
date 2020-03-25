for /f "tokens=1,2 delims=." %%a in ("%1") do (
	set ROOT_FILE=%%a
	set CII_EXTENTION=%%b
)
set NG=%2
set NX=%3
set NY=%4
set NZ=%5
set EQ=%6
set V1=%7
set V2=%8
cp -f source.dat Source_To_LU17/
cd Source_To_LU17/
main.exe %NG% %NX% %NY% %NZ% > source.txt
cd ..
cp -f Source_To_LU17/fort.17 %ROOT_FILE%.src
pre_cit %ROOT_FILE%.cii
mv %ROOT_FILE%.XSU XSU_MOD
cd XSU_MOD
neutyp %ROOT_FILE%.XSU %EQ% %V1% %V2% > xsu_data.dat
cd .. 
mv XSU_MOD/%ROOT_FILE%.XSU .
citvap0.exe %ROOT_FILE%.cii
caremdb -opt:export -val:meshflux %ROOT_FILE%.cdb
goto :EOF
