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

copy /Y /V x64\Profile\*.pgc x64\Release
@if %errorlevel% neq 0 exit /b %errorlevel%
copy /Y /V x64\Profile\*.pgd x64\Release
@if %errorlevel% neq 0 exit /b %errorlevel%

:: Build release exe

msbuild "Prime Bitstrings.vcxproj" -p:Configuration=Release -p:Platform=x64 -t:Build
