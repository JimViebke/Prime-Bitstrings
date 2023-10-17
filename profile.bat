:: Remove existing profile data

@del default.profdata
@del default.profraw

:: Remove any exes from the Profile directory

@del /s x64\Profile\*.exe

:: Build profile exe

msbuild "Prime Bitstrings.vcxproj" -p:Configuration=Profile -p:Platform=x64 -t:Build
@if %errorlevel% neq 0 exit /b %errorlevel%

:: Rename binary

@rename "x64\Profile\Prime Bitstrings.exe" "PB profile.exe"
@if %errorlevel% neq 0 exit /b %errorlevel%

:: Run profile exe

@start /min "" "x64\Profile\PB profile.exe"
@if %errorlevel% neq 0 exit /b %errorlevel%

@for /L %%i in (1, 1, 100) do @(	
	@tasklist | find /i "PB profile.exe" > NUL && echo Profiling... || goto next
	@timeout /t 1 /nobreak > NUL
)

@echo Profiling time limit reached
@taskkill /fi "IMAGENAME eq PB profile.exe" > NUL
@if %errorlevel% neq 0 exit /b %errorlevel%

:next
@echo Finished profiling

:: Convert profile data

llvm-profdata merge -output=default.profdata default.profraw
@if %errorlevel% neq 0 exit /b %errorlevel%

:: Build release exe

msbuild "Prime Bitstrings.vcxproj" -p:Configuration=Release -p:Platform=x64 -t:Build
@if %errorlevel% neq 0 exit /b %errorlevel%

:: Copy "current" -> "previous"
:: Copy new build -> "current"

@echo Copying 'current_profiled' build over 'previous_profiled' build
@copy /B /V /Y x64\Release\current_profiled.exe x64\Release\previous_profiled.exe

@echo Copying new build over 'current_profiled' build
@copy /B /V /Y x64\Release\"Prime Bitstrings.exe" x64\Release\current_profiled.exe
