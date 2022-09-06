
rem 
echo %*

echo n1=%~nx1
echo x1=%~x1

echo 1=%1
echo 2=%2
echo 3=%3

rem pause

j:\consalpha 80
rem cd /D %1\..\

rem call condabin\conda  activate  %2
r @call conda activate %2
cls

@echo anaconda environment: %2
@echo run=%3


@set pyname="%~nx1"

@rem @set activate=%~dp1\Scripts\activate

@if   "%ltx_conda%" == ""   set ltx_conda=conda

@rem @call "%activate%" %2
@call %ltx_conda% activate  %2
if %errorlevel% NEQ 0 (
  exit 172
)


cd
@echo on

@rem @python %3
@rem echo %pyname% %3
@rem cmd.exe

@if  "%debug.python.hook%" == "" (
            @%pyname% %3
    ) else (	
       @echo debug.python.hook=%debug.python.hook%
       @"%debug.python.hook%" %pyname% %3
       exit
    )
