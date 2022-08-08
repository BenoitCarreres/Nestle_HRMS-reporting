C:\R\R-4.0.2\bin\Rscript.exe allData2report.R %*
IF ERRORLEVEL 1 (
  echo **************************  Reporting FAILURE  **************************
  PAUSE
)
