import geatpy
import pytest
import numpy as np


@geatpy.Problem.single
def eval_function(indiv):
    return np.sum(indiv)


@geatpy.Problem.single
def eval_function_with_cv(indiv):
    return np.sum(indiv), np.asarray(indiv)


def test_Problem_single():
    indivs = [[1, 2, 3], [4, 5, 6]]
    assert np.array_equal(eval_function(indivs), [[6], [15]])
    obj, cv = eval_function_with_cv(indivs)
    assert np.array_equal(obj, [[6], [15]])
    assert np.array_equal(cv, [[1, 2, 3], [4, 5, 6]])


@pytest.fixture
def population_only_with_Phen():
    indivs = np.array([[1, 2, 3], [4, 5, 6]])
    yield geatpy.Population(None, Phen=indivs)


def test_Problem_evaluation_aimFunc_evalVars_invalid(
        population_only_with_Phen):
    pop = population_only_with_Phen
    problem = geatpy.Problem('test',
                             M=1,
                             maxormins=[1],
                             Dim=3,
                             varTypes=np.zeros(3),
                             lb=np.zeros(3),
                             ub=10 * np.ones(3))
    with pytest.raises(RuntimeError):
        problem.evaluation(pop)


def test_Problem_evaluation_with_evalVars(population_only_with_Phen):
    pop = population_only_with_Phen
    problem = geatpy.Problem('test',
                             M=1,
                             maxormins=[1],
                             Dim=3,
                             varTypes=np.zeros(3),
                             lb=np.zeros(3),
                             ub=10 * np.ones(3),
                             evalVars=eval_function)
    problem.evaluation(pop)
    assert np.array_equal(pop.ObjV, [[6], [15]])


def test_Problem_evaluation_with_aimFunc_and_evalVars(
        population_only_with_Phen):

    def aim_func_example(pop):
        pop.ObjV = [[10]]

    pop = population_only_with_Phen
    problem = geatpy.Problem('test',
                             M=1,
                             maxormins=[1],
                             Dim=3,
                             varTypes=np.zeros(3),
                             lb=np.zeros(3),
                             ub=10 * np.ones(3),
                             aimFunc=aim_func_example,
                             evalVars=eval_function)
    problem.evaluation(pop)
    assert np.array_equal(pop.ObjV, [[10]])
