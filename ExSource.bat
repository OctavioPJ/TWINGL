for /f "tokens=1,2 delims=." %%a in ("%1") do (
	set ROOT_FILE=%%a
	set CII_EXTENTION=%%b
)
set NG=%2
echo %ROOT_FILE%
cp -f source.dat Source_To_LU17/
cd Source_To_LU17
main.exe %NG% > source.txt
cd ..
cp -f Source_To_LU17/fort.17 %ROOT_FILE%.src
pre_cit %ROOT_FILE%.cii
python Map.py %ROOT_FILE%.cii
cp %ROOT_FILE%.bxsu XSU_MOD/
mv %ROOT_FILE%.XSU XSU_MOD/
cd XSU_MOD
Bmain.exe %ROOT_FILE%.XSU > xsu_aux.dat
cd ..
mv XSU_MOD/%ROOT_FILE%.XSU .
citvap0.exe %ROOT_FILE%.cii
caremdb -opt:export -val:meshflux %ROOT_FILE%.cdb
goto :EOF