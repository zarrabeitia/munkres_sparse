"""
Implementation of the sparse kuhn-munkres algorithm.
For a general description of the full munkres algorithm, see
http://csclab.murraystate.edu/bob.pilgrim/445/munkres.html
This implementation assumes a non-complete graph and stores
the elements in a sparse matrix, following the ideas outlined
in http://dl.acm.org/citation.cfm?id=1388969.1389035
(Sailor Assignment Problem)
"""

__author__ = "Luis A. Zarrabeitia"
__copyright__ = "Copyright 2011, Luis A. Zarrabeitia"
__credits__ = ["Luis A. Zarrabeitia"]

__license__ = "GPLv2 or later"
__version__ = "0.1"
__maintainer__ = "Luis. A. Zarrabeitia"
__email__ = "zarrabeitia at gmail dot com"
__status__ = "Production"
__url__ = "https://github.com/zarrabeitia/munkres_sparse"


__all__ = "munkres"

# Define the float value "inf". This is more architechture-independant than
# using float("inf")
INFINITY = 1e1000


class MunkresMatrix(object):
    """
    Auxiliary class.
    Stores the sparse representation of the munkres matrix.
    The rows and columns are remapped (to avoid having emtpy rows and columns).
    """
    def __init__(self, values):
        assert all(value >= 0 for i, j, value in values)
        bigvalue = 1e10
        rowindices = list(set(i for (i, j, value) in values))
        colindices = list(set(j for (i, j, value) in values))
        rowmap = dict((k, v) for v, k in enumerate(rowindices))
        colmap = dict((k, v) for v, k in enumerate(colindices))
        self.transposed = False
        self.values = [(rowmap[i], colmap[j], value)
                        for (i, j, value) in values]
        self.values.sort()
        self.rowmap = rowindices
        self.colmap = colindices
        self.nrows = len(rowindices)
        self.real_columns = len(colindices)
        self.ncols = len(colindices) + self.nrows
        # Ensure there is a feasible but very undesireable solution
        # covering the rows
        for i in xrange(self.nrows):
            self.values.append((i, self.real_columns + i, bigvalue))
        self.K = self.nrows
        self.rowindices = xrange(self.nrows)
        self.colindices = xrange(self.ncols)
        self.row_adds = [0] * self.nrows
        self.column_adds = [0] * self.ncols

    def remap(self, indices):
        """
        Transform the list of indices back to the original
        domain.
        """
        return [(self.rowmap[i], self.colmap[j])
                for i, j in indices if j < self.real_columns]

    def row(self, rowindex):
        """
        Returns the list of (value, column) in row rowindex
        """
        return ((value + self.column_adds[j] + self.row_adds[i], j)
                for i, j, value in self.values if i == rowindex)

    def get_values(self):
        """
        Returns the current values of the matrix.
        """
        return ((i, j, value + self.column_adds[j] + self.row_adds[i])
                for i, j, value in self.values)

    def add_column(self, colindex, value):
        """
        Adds value to all the elements of column colindex.
        """
        self.column_adds[colindex] += value

    def add_row(self, rowindex, value):
        """
        Adds value to all the elements of row rowindex.
        """
        self.row_adds[rowindex] += value

    def zeros(self):
        """
        Returns the indices (row, col) of all zero elements in the
        matrix. An element is considered to be zero if abs(value) <= 1e-6
        """
        return [(i, j)
                for (i, j, value) in self.get_values() if abs(value) <= 1e-6]


class Munkres(object):
    """
    Auxiliary class. Use the top level munkres method instead.
    """
    def __init__(self, values):
        """
        Initialize the munkres.
        values: list of non-infinite values entries of the cost matrix
                [(i,j,value)...]
        """
        self.matrix = MunkresMatrix(values)
        self.starred = set()
        self.primed = set()
        self.covered_columns = [False] * self.matrix.ncols
        self.covered_rows = [False] * self.matrix.nrows
        self.last_primed = None

    def munkres(self):
        """
        Executes the munkres algorithm.
        Returns the optimal matching.
        """
        next_step = self._step_1
        while next_step:
            next_step = next_step()

        # Transform the mapping back to the input domain
        return self.matrix.remap(self.starred)

    def _step_1(self):
        """
        For each row of the matrix, find the smallest element and subtract it
        from every element in its row.  Go to Step 2.
        """
        # TODO: This can probably be done much better than using .row(i),
        # but it is executed only once, so the performance penalty is low.
        for i in self.matrix.rowindices:
            minimum = min(self.matrix.row(i))[0]
            self.matrix.add_row(i, -minimum)
        return self._step_2

    def _step_2(self):
        """
        Find a zero (Z) in the resulting matrix.  If there is no starred zero
        in its row or column, star Z.
        Repeat for each element in the matrix. Go to Step 3.
        """
        zeros = self.matrix.zeros()
        for (i, j) in zeros:
            for (i1, j1) in self.starred:
                if (i1 == i or j1 == j):
                    break
            else:
                self.starred.add((i, j))
        return self._step_3

    def _step_3(self):
        """
        Cover each column containing a starred zero.  If K columns are covered,
        the starred zeros describe a complete set of unique assignments.  In
        this case, Go to DONE, otherwise, Go to Step 4.
        """
        for (_, j) in self.starred:
            self.covered_columns[j] = True
        if sum(self.covered_columns) == self.matrix.K:
            return None
        else:
            return self._step_4

    def _find_uncovered_zero(self):
        """
        Returns the (row, column) of one of the uncovered zeros in the matrix.
        If there are no uncovered zeros, returns None
        """
        zeros = self.matrix.zeros()
        for (i, j) in zeros:
            if not self.covered_columns[j] and not self.covered_rows[i]:
                return (i, j)
        return None

    def _step_4(self):
        """
        Find a noncovered zero and prime it.  If there is no starred zero in
        the row containing this primed zero, Go to Step 5.  Otherwise, cover
        this row and uncover the column containing the starred zero. Continue
        in this manner until there are no uncovered zeros left. Save the
        smallest uncovered value and Go to Step 6.
        """
        done = False
        while not done:
            zero = self._find_uncovered_zero()
            if zero:
                i, j = zero
                self.primed.add((i, j))
                self.last_primed = (i, j)
                st = [(i1, j1) for (i1, j1) in self.starred if i1 == i]
                if not st:
                    return self._step_5
                assert len(st) == 1
                i1, j1 = st[0]
                self.covered_rows[i] = True
                self.covered_columns[j1] = False
            else:
                done = True
        return self._step_6

    def _step_5(self):
        """
        Construct a series of alternating primed and starred zeros as follows.
        Let Z0 represent the uncovered primed zero found in Step 4. Let Z1
        denote the starred zero in the column of Z0 (if any). Let Z2 denote
        the primed zero in the row of Z1 (there will always be one). Continue
        until the series terminates at a primed zero that has no starred zero
        in its column.  Unstar each starred zero of the series, star each
        primed zero of the series, erase all primes and uncover every line in
        the matrix. Return to Step 3.
        """
        last_primed = self.last_primed
        last_starred = None
        primed = [last_primed]
        starred = []
        while True:
            # find the starred zero in the same column of last_primed
            t = [(i, j) for (i, j) in self.starred if j == last_primed[1]]
            if not t:
                break
            assert len(t) == 1
            last_starred = t[0]
            starred.append(last_starred)
            t = [(i, j) for (i, j) in self.primed if i == last_starred[0]]
            assert len(t) == 1
            last_primed = t[0]
            primed.append(last_primed)
        for s in starred:
            self.starred.remove(s)
        for p in primed:
            self.starred.add(p)
        self.primed.clear()
        for i in xrange(len(self.covered_rows)):
            self.covered_rows[i] = False

        return self._step_3

    def _step_6(self):
        """
        Add the value found in Step 4 to every element of each covered row, and
        subtract it from every element of each uncovered column.  Return to
        Step 4 without altering any stars, primes, or covered lines.
        """
        minval = INFINITY
        for i, j, value in self.matrix.get_values():
            covered = self.covered_rows[i] or self.covered_columns[j]
            if not covered and minval > value:
                minval = value
        assert 1e-6 < abs(minval) < INFINITY
        for i in self.matrix.rowindices:
            if self.covered_rows[i]:
                self.matrix.add_row(i, minval)
        for j in self.matrix.colindices:
            if not self.covered_columns[j]:
                self.matrix.add_column(j, -minval)
        return self._step_4

import random
import itertools


def random_test_munkres(nrows, ncols):
    """
    Naive test for the munkres implementation.
    Generates a random sparse cost matrix, applies munkres, and compares the
    result with the exahustive search.
    nrows, ncols: number of rows and columns of the generated matrix
    """
    values = [(i, j, random.random())
                for i in xrange(nrows)
                for j in xrange(ncols) if random.random() > .8]
    values_dict = dict(((i, j), v) for i, j, v in values)
    print values
    munkres_match = munkres(values)
    munkres_weight = sum(values_dict[p] for p in munkres_match)
    print len(munkres_match)
    print munkres_match
    print munkres_weight
    minimum = min(nrows, ncols)
    rows = set(i for i, j, v in values)
    cols = set(j for i, j, v in values)
    for part_row in itertools.combinations(rows, minimum):
        for part_col in itertools.combinations(cols, minimum):
            matching = zip(part_row, part_col)
            weight = sum(values_dict.get(p, INFINITY) for p in matching)
            if weight < munkres_weight:
                print "Munkres failed"
                print values
                print weight
                print matching
                print "munkres weight", munkres_weight
                raise Exception()
    return munkres_weight


def munkres(costs):
    """
    Entry method to solve the assignment problem.
    costs: list of non-infinite values entries of the cost matrix
            [(i,j,value)...]
    """
    solver = Munkres(costs)
    return solver.munkres()
