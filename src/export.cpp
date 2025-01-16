#include "export.hpp"

void exportResults(Grid *grid, const PiercedVector<double>& levelSetValues, std::string fileName, int iteration) 
{
    std::string file;
    file = fileName + std::to_string(iteration);
    grid->setExportName(file);

    std::vector<double> LSV = PVtoV(levelSetValues, grid);
    grid->addData("LevelSet", LSV);

    grid->write();
}
