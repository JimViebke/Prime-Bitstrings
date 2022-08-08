:: Remove any existing PGO profiles (.pgc) and profile databases (.pgd)

del /s x64\*.pgc
del /s x64\*.pgd

:: Clean both targets

:: msbuild "Prime Bitstrings.vcxproj" -p:Configuration=Profile -p:Platform=x64 -t:Clean
:: @if %errorlevel% neq 0 exit /b %errorlevel%
:: msbuild "Prime Bitstrings.vcxproj" -p:Configuration=Release -p:Platform=x64 -t:Clean
:: @if %errorlevel% neq 0 exit /b %errorlevel%

:: Build profile exe

msbuild "Prime Bitstrings.vcxproj" -p:Configuration=Profile -p:Platform=x64 -t:Build
@if %errorlevel% neq 0 exit /b %errorlevel%

:: Run profile exe

x64\Profile\"Prime Bitstrings.exe"
@if %errorlevel% neq 0 exit /b %errorlevel%

:: Copy profile results

copy /V /Y x64\Profile\*.pgc x64\Release
@if %errorlevel% neq 0 exit /b %errorlevel%
copy /V /Y x64\Profile\*.pgd x64\Release
@if %errorlevel% neq 0 exit /b %errorlevel%

:: Build release exe

msbuild "Prime Bitstrings.vcxproj" -p:Configuration=Release -p:Platform=x64 -t:Build
@if %errorlevel% neq 0 exit /b %errorlevel%

:: Copy "current" -> "previous"
:: Copy new build -> "current"

@echo Copying 'current_profiled' build over 'previous_profiled' build
@copy /B /V /Y x64\Release\current_profiled.exe x64\Release\previous_profiled.exe
@if %errorlevel% neq 0 exit /b %errorlevel%

@echo Copying new build over 'current_profiled' build
@copy /B /V /Y x64\Release\"Prime Bitstrings.exe" x64\Release\current_profiled.exe
@if %errorlevel% neq 0 exit /b %errorlevel%
