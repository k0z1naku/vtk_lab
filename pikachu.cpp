#include <iostream>
#include <cmath>
#include <vector>
#include <cassert>
#include <algorithm>

#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkTetra.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkSmartPointer.h>

#include <gmsh.h>

using namespace std;

const double EAR_MOVEMENT_AMPLITUDE = 0.9;
const double TAIL_MOVEMENT_AMPLITUDE = 0.7; 

class CalcNode
{
friend class CalcMesh;

protected:
    double x;
    double y;
    double z;
    double smth;
    double vx;
    double vy;
    double vz;

public:
    CalcNode() : x(0.0), y(0.0), z(0.0), smth(0.0), vx(0.0), vy(0.0), vz(0.0) {}

    CalcNode(double x, double y, double z, double smth, double vx, double vy, double vz) 
        : x(x), y(y), z(z), smth(smth), vx(vx), vy(vy), vz(vz) {}

    void move(double tau) {
        x += vx * tau;
        y += vy * tau;
        z += vz * tau;
    }
};

class Element
{
friend class CalcMesh;

protected:
    unsigned long nodesIds[4];
};

class CalcMesh
{
protected:
    vector<CalcNode> nodes;
    vector<Element> elements;
    vector<int> earNodeIndices; 
    vector<int> tailNodeIndices; 

public:
    CalcMesh(const std::vector<double>& nodesCoords, const std::vector<std::size_t>& tetrsPoints, 
             const std::vector<int>& earIndices, const std::vector<int>& tailIndices) 
        : earNodeIndices(earIndices), tailNodeIndices(tailIndices) {

        nodes.resize(nodesCoords.size() / 3);
        for(unsigned int i = 0; i < nodesCoords.size() / 3; i++) {
            double pointX = nodesCoords[i*3];
            double pointY = nodesCoords[i*3 + 1];
            double pointZ = nodesCoords[i*3 + 2];
            double smth = pow(pointX, 2) + pow(pointY, 2) + pow(pointZ, 2);
            nodes[i] = CalcNode(pointX, pointY, pointZ, smth, 0.0, 0.0, 0.0);
        }

        elements.resize(tetrsPoints.size() / 4);
        for(unsigned int i = 0; i < tetrsPoints.size() / 4; i++) {
            elements[i].nodesIds[0] = tetrsPoints[i*4] - 1;
            elements[i].nodesIds[1] = tetrsPoints[i*4 + 1] - 1;
            elements[i].nodesIds[2] = tetrsPoints[i*4 + 2] - 1;
            elements[i].nodesIds[3] = tetrsPoints[i*4 + 3] - 1;
        }
    }

    void doTimeStep(double tau, unsigned int step) {
        for(unsigned int i = 0; i < nodes.size(); i++) {
            // Проверяем, является ли текущий узел одним из узлов ушей
            if (std::find(earNodeIndices.begin(), earNodeIndices.end(), i) != earNodeIndices.end()) {
                // Двигаем уши с учетом амплитуды
                nodes[i].x += sin(step * tau) * EAR_MOVEMENT_AMPLITUDE;
                nodes[i].y += cos(step * tau) * EAR_MOVEMENT_AMPLITUDE;
            }
            // Проверяем, является ли текущий узел одним из узлов хвоста
            else if (std::find(tailNodeIndices.begin(), tailNodeIndices.end(), i) != tailNodeIndices.end()) {
                // Двигаем хвост с учетом амплитуды
                nodes[i].x += cos(step * tau) * TAIL_MOVEMENT_AMPLITUDE;
                nodes[i].z += sin(step * tau) * TAIL_MOVEMENT_AMPLITUDE;
            }
            // В остальных случаях двигаем обычным образом
            else {
                nodes[i].move(tau);
            }
        }
    }

    void snapshot(unsigned int snap_number) {
        vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
        vtkSmartPointer<vtkPoints> dumpPoints = vtkSmartPointer<vtkPoints>::New();
        auto smth = vtkSmartPointer<vtkDoubleArray>::New();
        smth->SetName("smth");
        auto vel = vtkSmartPointer<vtkDoubleArray>::New();
        vel->SetName("velocity");
        vel->SetNumberOfComponents(3);

        for(unsigned int i = 0; i < nodes.size(); i++) {
            dumpPoints->InsertNextPoint(nodes[i].x, nodes[i].y, nodes[i].z);
            double _vel[3] = {nodes[i].vx, nodes[i].vy, nodes[i].vz};
            vel->InsertNextTuple(_vel);
            smth->InsertNextValue(nodes[i].smth);
        }

        unstructuredGrid->SetPoints(dumpPoints);
        unstructuredGrid->GetPointData()->AddArray(vel);
        unstructuredGrid->GetPointData()->AddArray(smth);

        for(unsigned int i = 0; i < elements.size(); i++) {
            auto tetra = vtkSmartPointer<vtkTetra>::New();
            tetra->GetPointIds()->SetId( 0, elements[i].nodesIds[0] );
            tetra->GetPointIds()->SetId( 1, elements[i].nodesIds[1] );
            tetra->GetPointIds()->SetId( 2, elements[i].nodesIds[2] );
            tetra->GetPointIds()->SetId( 3, elements[i].nodesIds[3] );   
            unstructuredGrid->InsertNextCell(tetra->GetCellType(), tetra->GetPointIds());
        }

        string fileName = "pikachu-step-" + std::to_string(snap_number) + ".vtu";
        vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
        writer->SetFileName(fileName.c_str());
        writer->SetInputData(unstructuredGrid);
        writer->Write();
    }
};

int main()
{
    double h = 4.0;
    double tau = 1;

    const unsigned int GMSH_TETR_CODE = 4;

    gmsh::initialize();
    gmsh::model::add("pikachu");

    try {
        gmsh::merge("pikachu.stl");
    } catch(...) {
        gmsh::logger::write("Could not load STL mesh: bye!");
        gmsh::finalize();
        return -1;
    }

    double angle = 40;
    bool forceParametrizablePatches = false;
    bool includeBoundary = true;
    double curveAngle = 180;
    gmsh::model::mesh::classifySurfaces(angle * M_PI / 180., includeBoundary, forceParametrizablePatches, curveAngle * M_PI / 180.);
    gmsh::model::mesh::createGeometry();

    std::vector<std::pair<int, int> > s;
    gmsh::model::getEntities(s, 2);
    std::vector<int> sl;
    for(auto surf : s) sl.push_back(surf.second);
    int l = gmsh::model::geo::addSurfaceLoop(sl);
    gmsh::model::geo::addVolume({l});

    gmsh::model::geo::synchronize();

    int f = gmsh::model::mesh::field::add("MathEval");
    gmsh::model::mesh::field::setString(f, "F", "4");
    gmsh::model::mesh::field::setAsBackgroundMesh(f);

    gmsh::model::mesh::generate(3);

    std::vector<double> nodesCoord;
    std::vector<std::size_t> nodeTags;
    std::vector<double> parametricCoord;
    gmsh::model::mesh::getNodes(nodeTags, nodesCoord, parametricCoord);

    std::vector<std::size_t>* tetrsNodesTags = nullptr;
    std::vector<int> elementTypes;
    std::vector<std::vector<std::size_t>> elementTags;
    std::vector<std::vector<std::size_t>> elementNodeTags;
    gmsh::model::mesh::getElements(elementTypes, elementTags, elementNodeTags);
    for(unsigned int i = 0; i < elementTypes.size(); i++) {
        if(elementTypes[i] != GMSH_TETR_CODE)
            continue;
        tetrsNodesTags = &elementNodeTags[i];
    }

    if(tetrsNodesTags == nullptr) {
        cout << "Can not find tetra data. Exiting." << endl;
        gmsh::finalize();
        return -2;
    }

    cout << "The model has " <<  nodeTags.size() << " nodes and " << tetrsNodesTags->size() / 4 << " tetrs." << endl;

    for(int i = 0; i < nodeTags.size(); ++i) {
        assert(i == nodeTags[i] - 1);
    }
    assert(tetrsNodesTags->size() % 4 == 0);

    std::vector<int> earIndices = { 23,  76}; 
    std::vector<int> tailIndices = {10, 26, 27, 28, 29, 31, 32, 33, 34, 35, 36, 37, 38, 39}; 

    CalcMesh mesh(nodesCoord, *tetrsNodesTags, earIndices, tailIndices);

    gmsh::finalize();

    mesh.snapshot(0);
    for(unsigned int step = 1; step < 100; step++) {
        mesh.doTimeStep(tau, step);
        mesh.snapshot(step);
    }

    return 0;
}
