@echo off

set a=x64\Release\previous_profiled.exe
set b=x64\Release\current_profiled.exe

%a% && %b% && %a% && %b% && %a% && %b%
