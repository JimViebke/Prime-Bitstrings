@echo off

set a=x64\Release\previous_profiled.exe benchmark
set b=x64\Release\current_profiled.exe benchmark

:loop

%a% && %b% && %a% && %b%

@goto loop
