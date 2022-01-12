#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "ContactInhibitionCellCycleModel_g.hpp"
#include "VolumeTrackingModifier.hpp"
#include "WildTypeCellMutationState.hpp"

#include "AbstractCellBasedSimulationModifier.hpp"
#include "AbstractForce.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "NodesOnlyMesh.hpp"
#include "CellsGenerator.hpp"
#include "UniformCellCycleModel.hpp"
#include "TransitCellProliferativeType.hpp"
#include "GeneralisedLinearSpringForce_g.hpp"
#include "OffLatticeSimulation.hpp"
#include "SmartPointers.hpp"
#include "DiffusionForce_g.hpp"
#include "RandomNumberGenerator.hpp"
//This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

#include <NodeAttributes.hpp>

#include "RandomCellKiller.hpp"

#include <iostream>
#include <random>
#include <math.h>

#include <AbstractTwoBodyInteractionForce.hpp>


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class CellCounter : public AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>
{
private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellWriter<ELEMENT_DIM, SPACE_DIM> >(*this);
    }
	unsigned int count = 0;
public:

    CellCounter()
        : AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>("CellCounts.dat")
    {
    }

    void VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
    {

	if (fmod(SimulationTime::Instance()->GetTime(), .25) == 0)
	{
    c_vector<double, SPACE_DIM> cell_location = pCellPopulation->GetLocationOfCellCentre(pCell);
		count++;

    for (unsigned i=0; i<SPACE_DIM; i++)
    {
    *this->mpOutStream << cell_location[i] << " ";
    }

		if(count == pCellPopulation->GetNumAllCells())
		{
      // *this->mpOutStream << "\n";

			// *this->mpOutStream << pCellPopulation->GetNumAllCells() << "\n";

			count = 0;
		}
	}

   }
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(CellCounter)
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(CellCounter)

class TestTumorSimulation : public AbstractCellBasedTestSuite
{
  public:

void TestTumor()
{
  HoneycombMeshGenerator generator(20, 20);
  MutableMesh<2,2>* p_generating_mesh = generator.GetCircularMesh(10);
  NodesOnlyMesh<2> mesh;
  mesh.ConstructNodesWithoutMesh(*p_generating_mesh, 1.5);

  std::vector<CellPtr> cells;
  MAKE_PTR(TransitCellProliferativeType, p_transit_type);
  MAKE_PTR(WildTypeCellMutationState, p_state);

  RandomNumberGenerator::Instance()->Reseed(300);

  std::normal_distribution<double> G1(10.7, 3.0);
  std::normal_distribution<double> S(6.4, 3.0);
  std::normal_distribution<double> G2(4.2, .4);
  std::normal_distribution<double> M(1.0, .1);
  std::default_random_engine gen(300);

	for (unsigned i = 0; i<p_generating_mesh->GetNumNodes(); i++)
  // for (unsigned i = 0; i<1; i++)
	{
		ContactInhibitionCellCycleModel* p_cycle_model = new ContactInhibitionCellCycleModel();
		p_cycle_model->SetDimension(2);
		p_cycle_model->SetBirthTime((RandomNumberGenerator::Instance()->ranf() * 40) - 25);
		p_cycle_model->SetQuiescentVolumeFraction(0.5);
		p_cycle_model->SetEquilibriumVolume(1.0);

    p_cycle_model->SetTransitCellG1Duration(G1(gen));
    p_cycle_model->SetSDuration(S(gen));
    p_cycle_model->SetG2Duration(G2(gen));
    p_cycle_model->SetMDuration(M(gen));
		// p_cycle_model->SetTransitCellG1Duration(-10.7 * log(1.0001 - (RandomNumberGenerator::Instance()->ranf())));
		// p_cycle_model->SetSDuration(-6.4 * log(1.0001 - (RandomNumberGenerator::Instance()->ranf())));
		// p_cycle_model->SetG2Duration(-4.2 * log(1.0001 - (RandomNumberGenerator::Instance()->ranf())));
		// p_cycle_model->SetMDuration(-1 * log(1.0001 - (RandomNumberGenerator::Instance()->ranf())));


		CellPtr p_cell(new Cell(p_state, p_cycle_model));
		p_cell->SetCellProliferativeType(p_transit_type);
		p_cell->InitialiseCellCycleModel();

		cells.push_back(p_cell);
	}

  NodeBasedCellPopulation<2> cell_population(mesh, cells);

  cell_population.AddCellWriter<CellCounter>();

  OffLatticeSimulation<2> simulator(cell_population);
  simulator.SetOutputDirectory("Python_Chaste");
  simulator.SetSamplingTimestepMultiple(10);
  simulator.SetEndTime(75);

  MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
  simulator.AddForce(p_force);
  p_force->SetCutOffLength(1.5);

  MAKE_PTR(DiffusionForce<2>, d_force);
  d_force->SetAbsoluteTemperature(.1);
  simulator.AddForce(d_force);

  MAKE_PTR(VolumeTrackingModifier<2>, v_modifier);
  simulator.AddSimulationModifier(v_modifier);

  MAKE_PTR_ARGS(RandomCellKiller<2>, p_cell_killer, (&cell_population, 0.007));
  simulator.AddCellKiller(p_cell_killer);

  simulator.Solve();
  }
};



/* Code for copying into Chaste
docker cp ~/Desktop/Chaste_Code/graydon.hpp 04f64df5bf1a:/home/chaste/projects/graydon/test
docker cp ~/Desktop/Chaste_Code/DiffusionForce.cpp 04f64df5bf1a:/home/chaste/src/cell_based/src/population/forces/
docker cp ~/Desktop/Chaste_Code/ContactInhibitionCellCycleModel.cpp 04f64df5bf1a:/home/chaste/src/cell_based/src/cell/cycle/
docker cp ~/Desktop/Chaste_Code/AbstractSimplePhaseBasedCellCycleModel.cpp 04f64df5bf1a:/home/chaste/src/cell_based/src/cell/cycle/


docker cp ~/Desktop/src 04f64df5bf1a:/home/chaste/projects/graydon/src



docker cp ~/Desktop/Chaste_Code/GeneralisedLinearSpringForce.cpp 04f64df5bf1a:/home/chaste/src/cell_based/src/population/forces/

docker cp 04f64df5bf1a:/home/chaste/testoutput ~/Desktop/Chaste/

docker cp 04f64df5bf1a:/home/chaste/src/cell_based/src/population/forces/GeneralisedLinearSpringForce.cpp ~/Desktop/Chaste_Code/



*/
