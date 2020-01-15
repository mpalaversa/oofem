# this example illustrates an injection of external temperature field during runtime
import oofempy
import util

#Function returning temperature.
def evalField(coords, mode, tStep):
    val = 10.*tStep.giveNumber()*coords[0]
    print("Evaluating field at %f,%f. Assigning temperature %f" % (coords[0], coords[1], val))
    return (val,) #Return list of len 1


def test_5():
    # engngModel
    problem = oofempy.staticStructural(nSteps=3, outFile="test_5.out")

    # domain (if no engngModel specified to domain, it is asigned to the last one created)
    domain = oofempy.domain(1, 1, problem, oofempy.domainType._2dPlaneStressMode, tstep_all=True, dofman_all=True, element_all=True)
    problem.setDomain(1, domain, True)
    
    ltf1 = oofempy.constantFunction(1, domain, f_t=1.)
    ltf2 = oofempy.piecewiseLinFunction(2, domain, t=(1., 5.), f_t= (0., 4.))
    ltfs = (ltf1,ltf2)

    # boundary conditions
    bc1 = oofempy.boundaryCondition(1, domain, loadTimeFunction=1, prescribedValue=0.0)
    n2 = oofempy.nodalLoad(2, domain, loadTimeFunction=2, components=(0.,0.), dofs=(1,2))
    bcs = (bc1, n2)

    # nodes
    n1 = oofempy.node(1, domain, coords=(0.,  0., 0. ), bc=(bc1,bc1))
    n2 = oofempy.node(2, domain, coords=(2.4, 0., 0. ), bc=(0,bc1))
    n3 = oofempy.node(3, domain, coords=(2.4, 0.8, 0.), load=(n2,) )
    nodes = (n1,n2,n3)

    # material and cross section
    mat = oofempy.isoLE(1, domain, d=1., E=30.e3, n=0.2, tAlpha=1.2e-5)
    cs  = oofempy.simpleCS(1, domain, thick=0.5)
    
    # elements
    e1 = oofempy.trPlaneStress2d(1, domain, nodes=(1,2,3), mat=1, crossSect=1)
    elems = (e1,)

    # setup domain
    util.setupDomain(domain, nodes, elems, (mat,), (cs,), bcs, ltfs, ())

    print("\nSolving problem")
    problem.checkProblemConsistency()
    problem.init()
    problem.postInitialize()
    problem.giveTimer().startTimer(oofempy.EngngModelTimerType.EMTT_AnalysisTimer)
    activeMStep = problem.giveMetaStep(1)
    problem.initMetaStepAttributes(activeMStep);
    
    #Create dummy temperature field, define from where to obtain the values
    f = oofempy.PythonField()
    f.setModuleName('test_5')
    f.setFunctionName('evalField')
    context = problem.giveContext()
    field_man = context.giveFieldManager()
    #register field so OOFEM knows of external temperature. In every iteration, it calls PythonField::evaluateAt() which propagates further to evalField()
    field_man.registerField(f, oofempy.FieldType.FT_Temperature)
    
    for timeStep in range(3):
        problem.preInitializeNextStep()
        problem.giveNextStep()
        currentStep = problem.giveCurrentStep()
        problem.initializeYourself( currentStep )
        problem.solveYourselfAt( currentStep )
        problem.updateYourself( currentStep )
        problem.terminate( currentStep )
        print("TimeStep %d finished" % (timeStep))
    problem.terminateAnalysis()

    ##check solution
    v3 = problem.giveUnknownComponent(oofempy.ValueModeType.VM_Total, problem.giveCurrentStep(False), domain, domain.giveDofManager(3).giveDofWithID(oofempy.DofIDItem.D_v))
    assert (round (v3-4.608e-4, 8) == 0), "Node 3 dof 2 displacement check failed"

    problem.terminateAnalysis()
    print("\nProblem solved")


if __name__ == "__main__":
    test_5()
