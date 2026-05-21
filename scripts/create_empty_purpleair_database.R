# Create an empty purpleair duckdb database table

library(dplyr)

# table structure
pa_table <- tibble(sensor_index = numeric(0),
                      last_seen = lubridate::POSIXct(),
                      confidence = numeric(0),
                      pm2.5 = numeric(0))

# Path and filename of database file to be created
dbdir <- "./purpleair_duckdb"

# database connection - this creates the database if it doesn't already exist
con <- DBI::dbConnect(duckdb::duckdb(), dbdir = dbdir)

# Create the empty table
DBI::dbCreateTable(con, "pa_realtime", pa_table)

### Optional

# Test that some data can be added 
k <- "{YOUR PA API KEY}" # preferably loaded from an environment variable
pa_acquire_realtime(dbdir = dbdir, key = k)

# Should get about 10,000 records with the default

# take a look at some data
glimpse(tbl(con, "pa_realtime"))
  
# close the connection
DBI::dbDisconnect(con)

