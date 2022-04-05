# Known Issues in the RSVG branch

## Differences Between Local and Production versions

When deploying to Cloud Run, some unknown issue causes the production and local
versions to function differently. The following differences have been noticed
between local and production versions:

1. The download genomes feature does not function as intended

## Compare Locations/Mutations Plot does not display properly

When sequences before 1960 are included in the Compare Locations Plot, the
mutation data does not display properly. This is because the logic is based
on Unix time and does not account for negative values. Even one negative value
in collection_week results in the plot breaking
