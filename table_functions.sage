"""
This file contains classes for computing isogeny classes of Abelian varieties over finite fields.

The computation rests on the Honda-Tate theorem, which gives a bijection between isogeny classes
of abelian varieties of dimension g over F_q and Weil q-polynomials of degree 2g satisfying certain
conditions.  Recall that f is a Weil q-polynomial if all of the complex roots of f have absolute value
equal to the square root of q.  For more details on the conditions, see the `simplepow_brauer_data`
attribute.

The computation is organized into stages.  Each stage adds a number of attributes, which are then
loaded by later stages.

  - Enumerate simple isogeny classes using the `WeilPolynomials` iterator over Weil q-polynomials
  - Construct all isogeny classes as products of simple ones
  - Compute base changes from F_q to F_{q^r}, filling in primitive models, twists and endomorphism algebras
  - Combine data from external sources (number of isomorphism classes, number of Jacobians, etc)
  - Collate into a single file for uploading to Postgres

Each stage is represented by a `Stage` instance, which produces a list of `Tasks` (usually one
for each g and q).  A task can depend on previous tasks being completed, and there is a
`Controller` object that starts tasks running once their prerequisites have been met.
Tasks usually write their output to a file, and create an additional file to signal completion.
The sequence of stages to be executed is controlled by a configuration file `config.ini`.

The actual computation of the data associated to a Weil q-polynomial is done by an `IsogenyClass`
instance.  Attributes of an isogeny class that correspond to columns in the database
are decorated with a type, which specifies how they are loaded from and saved to a file.
The file format used is precisely that required by the `copy_from` and `reload` methods
of the `PostgresTable` class in the LMFDB.  The output of this script is thus easily loaded
into the LMFDB.

AUTHORS:

  - (2016-05-11) Taylor Dupuy, Kiran S. Kedlaya, David Roe, Christelle Vincent
  - (2019-02-02) David Roe, Everett Howe, added more tests for principal polarizations and Jacobians
  - (2019-05-01) Major refactoring and documentation

EXAMPLES:

Currently, the lmfdb can only be imported from the lmfdb folder, and then this file must be attached
from the root-unitary folder.  The actual computation must be done from the lmfdb folder,
since new instances of the `PostgresDatabase` object must be created.
The following example assumes that the lmfdb and root unitary repositories are installed side-by-side::

    sage: cd lmfdb
    sage: from lmfdb import db
    sage: cd ../root-unitary
    sage: %attach table_functions.sage
    sage: IC = IsogenyClasses()
    sage: cd ../lmfdb
    sage: IC.run_parallel(4)

The stages to be executed are controlled by the configuration file `config.ini`.
In addition to containing global configuration options, it specifies which columns should be stored by each stage.
An example configuration is included in this repository, and we describe the different sections here.

[extent]
This section controls the range of g and q to be computed.  For any (g, q) included, all
(g',q) with g' < g and all (g, q') with q'|q must also be included.
  - `g1`, `g2`, etc -- give a list of `q` to be included for this dimension.
[plan]
This section controls which stages are actually run.  By editing this section you can control
which stages are computed.
  - `stages` -- a list of stages, which will be run in order
[dirs]
This section describes where on the file system results are stored.
  - `base` -- a base folder
  - `subdirs` -- a list of subfolders to be created.  If the default values are kept,
     the subfolders play the following roles:
    - `basic` -- holds the results of the GenerateSimple and GenerateAll stages
    - `endalg` -- holds the results of the Basechange stage
    - `external` -- location for adding external data files for input to the Combine stage
    - `complete` -- holds the results of the Combine stage
    - `logs` -- holds logs for recording progress of computational stages
[logging]
This section contains configuration options for logging
  - `logfile` -- a template for where logs are stored (e.g. logs/{name}{data}.log)
  - `logfrequency` -- the frequency with which a status report will be posted to the log file
  - `logheader` -- a format string for a message to be printed when a task starts or finishes
[StageGenerateSimple], [StageGenerateAll], [StageBasechange],...
Each stage defined in this file has a corresponding section in the configuration file, specifying
input and output parameters
  - `in0`, `in1`, ... -- a list of input files, usually with formatting slots for g and q
  - `out0`, `out1`, ... -- a list of output files, usually with formatting slots for g and q
  - `data0`, `data1`, ... -- for each output file, a list of attributes to be saved to that file

Fields we want to populate with an example

label: "2.9.ab_d"
polynomial: [1,-1,3,-9,81]
angle_numbers (doubles): [0.23756..., 0.69210...]
number_field: "4.0.213413.1"
p-rank: 1
slopes: [0,1/2,1/2,1]
A_counts: [75, 7125]
C_counts: [9, 87]
known_jacobian (0,1,-1): 1
decomposition: ["9.2.-1._3"]
pricipally_polarizable (0,1,-1): 1
Brauer Invariants: inv_v( End(A_{FFbar_q})_{QQ} )=(v(\pi)/v(q))*[QQ(pi)_{v}: QQ(pi): v\vert p place of QQ(\pi)], these are stored as elements of QQ.
Primitive models: 

TODO: Add size to StageCombine, fix errors revealed by StageCheck
"""

######################################################################################################

from sage.databases.cremona import cremona_letter_code, class_to_int # for the make_label function
from sage.misc.lazy_attribute import lazy_attribute
from sage.all import Polynomial
import json, os, re, sys, time
import psutil
opj, ope = os.path.join, os.path.exists
from copy import copy
from collections import defaultdict, Counter
from multiprocessing import Process
from string import letters
from datetime import datetime
from itertools import combinations_with_replacement, imap, izip_longest
from ConfigParser import ConfigParser
from lmfdb.backend.database import PostgresDatabase
from psycopg2.sql import SQL, Identifier, Literal
from cypari2.handle_error import PariError

int_re = re.compile(r'-?\d+')

try:
    # Add the location of weil_polynomials.pyx to the load path
    sys.path.append(os.path.dirname(os.path.realpath(__file__)))
except NameError:
    pass
load("weil_polynomials.pyx")

######################################################################################################

# Some combinatorial utilities for enumerating non-simple isogeny classes from simple ones

class CombosWithReplacement(object):
    """
    We augment itertools' combinations_with_replacement with a length
    """
    def __init__(self, L, m):
        self.L = L
        self.m = m
    def __len__(self):
        return binomial(len(self.L)+self.m-1,self.m)
    def __iter__(self):
        return combinations_with_replacement(self.L, self.m)


# These functions are used for computing base changes

@cached_function
def symfunc(i, r):
    """
    Returns a symmetric function expressing the ith coefficient of the rth base change map.

    Namely, if f(x) = 1 + ... + a_n x^n and g(x) = 1 + ... + b_n x^n is the polynomial
    whose roots are the rth powers of the roots of f, then b_i = symfunc(i, r)(1, a1, ..., a_n)

    EXAMPLES:

        sage: symfunc(1,2)
        e[1, 1] - 2*e[2]
        sage: symfunc(2,2)
        e[2, 2] - 2*e[3, 1] + 2*e[4]
    """
    Sym = SymmetricFunctions(QQ)
    p = Sym.powersum()
    if i == 0:
        return p.one()
    e = Sym.elementary()
    return e(p(e[i]).map_support(lambda A: Partition([r*c for c in list(A)])))

@cached_function
def basechange_transform(g, r):
    """
    Returns a transform that takes in a `q`-Weil L-polynomials of degree `2g` (constant coefficient 1)
    and returns a `q^r`-Weil L-polynomial of degree `2g` whose roots are the `r`th powers of the roots
    of the input.
    """
    f = [symfunc(i, r) for i in range(g+1)]
    coeffs = [b.coefficients() for b in f]
    exps = [[{a: list(elem).count(a) for a in set(elem)} for elem in sorted(b.support()) if list(elem) and max(elem) <= 2*g] for b in f]
    def bc(Lpoly, q):
        # Assume that Lpoly has constant coefficient 1.
        R = Lpoly.parent()
        signed_coeffs = [(-1)^j * c for j, c in enumerate(Lpoly)]
        bc_coeffs = [1]
        for i in range(1, g+1):
            bc_coeffs.append((-1)^i*sum(c*prod(signed_coeffs[j]^e for j,e in D.iteritems()) for c, D in zip(coeffs[i], exps[i])))
        for i in range(1,g+1):
            # a_{g+i} = q^(ri) * a_{g-i}
            bc_coeffs.append(q^(r*i) * bc_coeffs[g-i])
        return R(bc_coeffs)
    return bc

def tensor_charpoly(f, g):
    r"""
    INPUT:

    - ``f`` -- the characteristic polynomial of a linear transformation
    - ``g`` -- the characteristic polynomial of a linear transformation

    OUTPUT:

    The characteristic polynomial of the tensor product of the linear transformations

    EXAMPLES::

        sage: x = PolynomialRing(ZZ,"x").gen();
        sage: tensor_charpoly((x - 3) * (x + 2),  (x - 7) * (x + 5))
        (x - 21) * (x - 10) * (x + 14) * (x + 15)
    """
    R.<y> = g.parent()[]
    A = f(y)
    B = R(g.homogenize(y))
    return B.resultant(A)

def base_change(Lpoly, r, algorithm=None, g = None, q = None, prec=53):
    """
    Returns a polynomial whose roots are the `r`th powers of the roots of the input polynomial.

    INPUT:

    - ``Lpoly`` -- a Weil `q`-polynomial of degree `2g`
    - ``r`` -- a positive integer
    - ``algorithm`` -- either "sym", "approx" or "linalg".
      - If "sym", an algorithm based on symmetric polynomials is used.
        This allows for caching the transformation and fast evaluation, but has a rapidly growing
        memory footprint when `r` gets larger than about 10.
      - If "approx", roots are found over a complex field, raised to the `r`th power
        and the base change polynomial is reconstructed.  For large `r` you must give an
        appropriate precision manually.
      - If "linalg", the base change is computed by evaluating at roots of unity times
        the generator of the polynomial ring.  This is slower than "sym" but does not have
        the memory issues, and no precision analysis is needed since the arithmetic is exact.
      If no algorithm is given, `r` is factored and "sym" and "linalg" used recursively, with
      "linalg" used for primes larger than 7

    Returns a transform that takes in a `q`-Weil L-polynomials of degree `2g` (constant coefficient 1)
    and returns a `q^r`-Weil L-polynomial of degree `2g` whose roots are the `r`th powers of the roots
    of the input.
    """
    if g is None:
        g = Lpoly.degree()
        assert g % 2 == 0
        g = g // 2
    if q is None:
        q = Lpoly.leading_coefficient().nth_root(g)
    if algorithm is None:
        # It's better to divide by larger numbers if possible, as long as we don't get so big that
        # computing the basechange_transform starts to cause memory issues
        for ell in range(7,1,-1):
            while r % ell == 0:
                Lpoly = base_change(Lpoly, ell, algorithm='sym', g=g, q=q)
                q = q^ell
                r = r//ell
        for ell, e in ZZ(r).factor():
            for i in range(e):
                Lpoly = base_change(Lpoly, ell, algorithm='linalg', g=g, q=q)
        return Lpoly
    elif algorithm == 'approx':
        C = ComplexField(prec)
        R = RealField(prec)
        LC = Lpoly.change_ring(C)
        x = LC.parent().gen()
        approx = prod((1 - x/alpha^r)^e for alpha, e in LC.roots())
        approx_coeffs = approx.list()
        acceptable_error = R(2)^-(prec//2)
        exact_coeffs = [c.real().round() for c in approx_coeffs]
        if max(abs(ap - ex) for ap, ex in zip(approx_coeffs, exact_coeffs)) > acceptable_error:
            raise RuntimeError
        return Lpoly.parent()(exact_coeffs)
    elif algorithm == 'linalg':
        # From https://github.com/edgarcosta/endomorphisms/blob/master/endomorphisms/OverFiniteField/utils.py
        K.<zeta> = CyclotomicField(r)
        R = Lpoly.parent()
        f = Lpoly.change_ring(K)
        T = f.parent().gen()
        poly = 1
        for i in range(r):
            poly *= f(zeta^i  * T)
        L = poly.list()
        newf = [None]*(1 + f.degree())
        for i, ci in enumerate(L):
            if i % r != 0:
                if ci != 0:
                    raise RuntimeError("i = %s, ci = %s" % (i, ci))
            else:
                newf[i/r] = ZZ(ci)
        return R(newf)
    else:
        return basechange_transform(g, r)(Lpoly, q)

class Cyclotomic(UniqueRepresentation):
    """
    An object for determining which polynomials are cyclotomic, and their order if they are.
    """
    # Convention: phi(m) = n
    def __init__(self):
        self.polynomials = {}
        self.mbound = 1
        self.nbound = 0
    def order(self, f):
        """
        If f is a cyclotomic polynomial Phi_m, returns m.
        Otherwise, returns 0.
        """
        if f.degree() >= self.nbound:
            self.compute(f.degree())
        return self.polynomials.get(f, 0)
    def get_mbound(self, nbound):
        """
        
        """
        # https://math.stackexchange.com/questions/265397/inversion-of-the-euler-totient-function
        steps = 0
        n = nbound
        while n != 1:
            n = euler_phi(n)
            steps += 1
        return 2 * 3^steps
    def compute(self, nbound):
        mbound = self.get_mbound(nbound)
        for m in range(self.mbound, mbound):
            self.polynomials[cyclotomic_polynomial(m)] = m
        self.nbound = nbound
        self.mbound = mbound

# These functions extend cremona codes to negative integers (prepending an a)

def signed_cremona_letter_code(m):
    if m >= 0:
        return cremona_letter_code(m)
    else:
        return 'a' + cremona_letter_code(-m)

def signed_class_to_int(code):
    if code == 'a':
        return ZZ(0)
    elif code.startswith('a'):
        return -ZZ(class_to_int(code[1:]))
    else:
        return ZZ(class_to_int(code))

# The following classes support nicely loading and saving
# to disk in a format readable by postgres

class PGType(lazy_attribute):
    check=None # bound for integer types
    def _load(self, x):
        """
        Wraps :meth:`load` for the appropriate handling of NULLs.
        """
        if x != r'\N':
            return self.load(x)
    def load(self, x):
        """
        Takes a string from a file and returns the appropriate
        Sage object.

        Should be inverse to :meth:`save`

        This default function can be overridden in subclasses.
        """
        if int_re.match(x):
            return ZZ(x)
        elif x.startswith('{'):
            return sage_eval(x.replace('{','[').replace('}',']'))
        elif x.startswith('['): # support jsonb
            return sage_eval(x)
        else:
            return x

    def _save(self, x):
        """
        Wraps :meth:`save` for the appropriate handling of NULLs.
        """
        if x is None:
            return r'\N'
        else:
            return self.save(x)

    def save(self, x, recursing=False):
        """
        Takes a Sage object stored in this attribute
        and returns a string appropriate to write to a file
        that postgres can load.

        Should be inverse to :meth:`load`

        This default function can be overridden in subclasses.
        """
        if isinstance(x, (list, tuple)):
            if self.pg_type == 'jsonb':
                return '[' + ','.join(PGType.save(self, a) for a in x) + ']'
            else:
                return '{' + ','.join(PGType.save(self, a) for a in x) + '}'
        else:
            if self.check and (x >= self.check or x <= -self.check):
                raise ValueError("Out of bounds")
            if isinstance(x, basestring):
                return '"' + x + '"'
            else:
                return str(x)

class pg_text(PGType):
    pg_type = 'text'
    def load(self, x):
        return x
    def save(self, x):
        return x
class pg_smallint(PGType):
    check = 2**15-1
    pg_type = 'smallint'
class pg_integer(PGType):
    check = 2**31-1
    pg_type = 'integer'
class pg_bigint(PGType):
    check = 2**63-1
    pg_type = 'bigint'
class pg_numeric(PGType):
    # Currently only used to store large integers, so no decimal handling needed
    pg_type = 'numeric'
class pg_smallint_list(PGType):
    check = 2**15-1
    pg_type = 'smallint[]'
class pg_integer_list(PGType):
    check = 2**31-1
    pg_type = 'integer[]'
class pg_bigint_list(PGType):
    check = 2**63-1
    pg_type = 'bigint[]'
class pg_float8_list(PGType):
    # We use double precision since that's what postgres returns as the type (and thus what database.py checks for)
    pg_type = 'double precision[]'
class pg_text_list(PGType):
    pg_type = 'text[]'
class pg_numeric_list(PGType):
    pg_type = 'numeric[]'
class pg_boolean(PGType):
    pg_type = 'boolean'
    def load(self, x):
        if x in ['t','1']:
            return True
        elif x in ['f','0']:
            return False
        else:
            raise RuntimeError
    def save(self, x):
        if x:
            return 't'
        else:
            return 'f'
class _rational_list(PGType):
    """
    A list of rational numbers (or nested such lists), stored as lists of strings.
    """
    def load(self, x):
        def recursive_QQ(y):
            if isinstance(y, basestring):
                return QQ(y)
            else:
                return map(recursive_QQ, y)
        x = PGType.load(self, x)
        return recursive_QQ(x)
    def save(self, x):
        def recursive_str(y):
            if isinstance(y, list):
                return [recursive_str(z) for z in y]
            else:
                return str(y)
        x = recursive_str(x)
        return PGType.save(self, x)
class _rational_mults(PGType):
    """
    A non-nested list of rational numbers, where multiplicities
    are emphasized by appending a letter, e.g. ["1/2A", "1/2B", "2/3A"]
    """
    def load(self, x):
        x = sage_eval(x.translate(None, letters).replace('{','[').replace('}',']')) # Remove letters
        return [QQ(y) for y in x]
    def save(self, x):
        cntr = Counter()
        def stringify(x):
            res = str(x) + cremona_letter_code(cntr[x]).upper()
            cntr[x] += 1
            return res
        return PGType.save(self, [stringify(y) for y in x])
class pg_rational_list(_rational_list):
    pg_type = 'text[]'
class pg_rational_mults(_rational_mults):
    pg_type = 'text[]'
class pg_jsonb(PGType):
    pg_type = 'jsonb'

class PGOld(PGType):
    """
    Columns from the current LMFDB schema.

    Even if the column name and type hasn't changed,
    it's still useful to use this class so that StageDBLoad works correctly.
    """
    def __init__(self, func=None, new_name=None, switch=(lambda x,db:x)):
        if func is None:
            self.new_name = new_name
            self.switch = switch
        else:
            if not hasattr(self, 'new_name'):
                self.new_name = None
            if not hasattr(self, 'switch'):
                self.switch = lambda x,db: x
            self(func)
    def __call__(self, func):
        PGType.__init__(self, func)
        if self.__name__.startswith('old_'):
            self.__name__ = self.__name__[4:]
        if self.new_name is None:
            self.new_name = self.__name__
        return self

class pg_text_old(pg_text, PGOld):
    pass
class pg_smallint_old(pg_smallint, PGOld):
    pass
class pg_integer_old(pg_integer, PGOld):
    pass
class pg_jsonb_old(pg_jsonb, PGOld): # More specific when using rationals??
    pass
class pg_boolean_old(pg_boolean, PGOld):
    pass
class pg_rational_list_old(_rational_list, PGOld):
    pg_type = 'jsonb'
class pg_rational_mults_old(_rational_mults, PGOld):
    pg_type = 'jsonb'

class Stage(object):
    def __init__(self, controller, input, output):
        self.controller = controller
        self.input = input
        self.output = output

    @lazy_attribute
    def tasks(self):
        return [self.Task(g, q, self) for g, q in self.controller.gq]

class GenericTask(object):
    def __init__(self, g, q, stage):
        self.g, self.q, self.stage = g, q, stage
        self.kwds = {'g':g, 'q':q}
        self.logheader = stage.controller.logheader.format(g=g, q=q, name=stage.shortname)
    @staticmethod
    def _done(filename):
        return filename[:-4] + '.done'
    def _logdata(self):
        return '_%s_%s' % (self.g, self.q)
    @lazy_attribute
    def _logfile(self):
        stage = self.stage
        controller = stage.controller
        return controller.logfile_template.format(name=stage.shortname, data=self._logdata())
    def ready(self):
        return all(ope(self._done(data[-1])) for data in self.input_data)
    def done(self):
        return all(ope(filename) for filename in self.donefiles)
    def log(self, s):
        s = str(s).strip()
        s += " (%s)" % (datetime.utcnow() - self.t0)
        s += '\n'
        with open(self._logfile, 'a') as F:
            F.write(s)
    @lazy_attribute
    def input_data(self):
        """
        Default behavior can be overridden in subclass
        """
        return [(None, filename.format(**self.kwds)) for filename in self.stage.input]
    @lazy_attribute
    def donefiles(self):
        return [self._done(output.format(**self.kwds)) for output, attributes in self.stage.output]
    def save(self, *sources):
        stage = self.stage
        controller = stage.controller
        for (outfile, attributes), source in zip(stage.output, sources):
            outfile = outfile.format(**self.kwds)
            controller.save(outfile, source, attributes, self.t0)
    def _run(self, db=None):
        if db is None:
            db = PostgresDatabase()
        self.t0 = datetime.utcnow()
        self.run(db)
        self.log("Finished")
    def enumerate(self, sequence, start=0):
        freq = self.stage.controller.logfrequency
        try:
            slen = '/%s' % len(sequence)
        except Exception:
            slen = ''
        for i, item in enumerate(sequence, start):
            if i and i%freq == 0:
                self.log('%s%s' % (i, slen))
            yield item

class Controller(object):
    def __init__(self, worker_count=1, config=None):
        if config is None:
            if os.path.exists('config.ini'):
                config = os.path.abspath('config.ini')
            else:
                raise ValueError("Must have config.ini in directory or specify location")
        self.config_file = config
        self.cfgp = cfgp = ConfigParser()
        cfgp.read(config)
        # Create subdirectories if they do not exist
        basedir = os.path.abspath(os.path.expanduser(cfgp.get('dirs', 'base')))
        if not ope(basedir):
            os.makedirs(basedir)
        subdirs = [sub.strip() for sub in cfgp.get('dirs', 'subdirs').split(',')] + ['logs']
        for subdir in subdirs:
            if not ope(opj(basedir, subdir)):
                os.mkdir(opj(basedir, subdir))
        self._extra_init()

        # Create Stages
        stages = []
        self.logfile_template = opj(basedir, cfgp.get('logging', 'logfile'))
        self.logfrequency = int(cfgp.get('logging', 'logfrequency'))
        self.logheader = cfgp.get('logging', 'logheader') + ' '
        plan = ['Stage'+stage.strip() for stage in cfgp.get('plan', 'stages').split(',')]
        for stage in plan:
            if stage not in cfgp.sections():
                raise ValueError("Missing section %s"%stage)
            info = [(key, cfgp.get(stage, key)) for key in cfgp.options(stage)]
            input = [opj(basedir, val) for (key, val) in info if key.startswith('in')]
            if any(key.startswith('out') for key, val in info):
                output, output_indexes = zip(*((opj(basedir, val), key[3:]) for (key, val) in info if key.startswith('out')))
            else:
                output_indexes = output = []
            if any(key.startswith('data') for key, val in info):
                data, data_indexes = zip(*(([attr.strip() for attr in val.split(',')], key[4:]) for (key, val) in info if key.startswith('data')))
            else:
                data_indexes = data = []
            if output_indexes != data_indexes:
                raise ValueError("Output and data specifications need to be in the same order for %s.\nOne was %s while the other was %s" % (stage, ', '.join(output_indexes), ', '.join(data_indexes)))
            output = zip(output, data)
            stages.append(getattr(self.__class__, stage)(self, input=input, output=output))

        self.stages = stages
        self.tasks = sum((stage.tasks for stage in stages), [])

    def _extra_init(self):
        # Parse config options before creating stages and tasks
        pass

    def clear_done(self):
        """
        .done files are used to incdicate that a task is complete
        (so that later tasks are able to start, and so that you can
        restart after a keyboard interrupt without losing all progress)

        Use this function to delete all such files and start from scratch.
        """
        for task in self.tasks:
            for filename in task.donefiles:
                if ope(filename):
                    os.unlink(filename)

    def run_serial(self, db=None):
        if db is None:
            db = PostgresDatabase()
        for task in self.tasks:
            if task.done():
                print "Already complete " + task.logheader
            else:
                print "Starting " + task.logheader
                task._run(db)

    def run_parallel(self, worker_count=8):
        processes = {}
        pcounter = 0
        tasks = copy(self.tasks)
        while tasks or processes:
            # start processes
            i = 0
            while i < len(tasks) and len(processes) < worker_count:
                if tasks[i].done():
                    task = tasks.pop(i)
                    print "Already complete " + task.logheader
                elif tasks[i].ready():
                    task = tasks.pop(i)
                    P = Process(target=task._run)
                    print "Starting " + task.logheader
                    P.start()
                    processes[pcounter] = (P, task)
                    pcounter += 1
                else:
                    i += 1
            time.sleep(0.1)
            # finish processes
            for p, (P, task) in list(processes.items()):
                if not P.is_alive():
                    processes.pop(p)
                    print "Finished " + task.logheader

class IsogenyClasses(Controller):
    """
    This class controls the creation of isogeny classes, grouped by g and q,
    as well as the process of loading and saving them to disk and to the LMFDB

    INPUT:

    - ``worker_count`` -- the number of processes to be allocated to this computation
    - ``config`` -- the filename for the configuration file, an example of which follows:

      [extent]
      g1 = 2-10,16
      g2 = 2-5
      [stage1]
      __dir__ = simple/
      all = all/
      complete = complete/

    Directories can be relative or absolute, and are indicated by ``__dir__``.
    Data specified within each stage is cumulative.
    """
    def __init__(self, worker_count=1, config=None):
        self.default_cls = IsogenyClass
        Controller.__init__(self, worker_count, config)


    def _extra_init(self):
        # Parse the extent section
        cfgp = self.cfgp
        gs = sorted([ZZ(gx[1:]) for gx in cfgp.options('extent')])
        if not gs:
            raise ValueError('Config file must specify [extent] section with lines like "g3=2-16,25"')
        self.gq = []
        gq_dict = {}
        for g in gs:
            qs_raw = cfgp.get('extent', 'g%s'%g)
            qs_split = qs_raw.split(',')
            qs = []
            for qrange in qs_split:
                qrange = qrange.strip()
                if qrange.count('-') == 1:
                    a,b = map(ZZ,qrange.split('-'))
                    for q in srange(a, b+1):
                        if q > 1 and q.is_prime_power():
                            qs.append(q)
                else:
                    q = ZZ(qrange)
                    if q.is_prime_power():
                        qs.append(q)
                    else:
                        raise ValueError("q=%s in g%s line is not a prime power" % (q, g))
            if not qs:
                raise ValueError("No qs specified for the given g")
            gq_dict[g] = D = set(qs)
            for q in qs:
                for d in q.divisors():
                    if d != 1 and d not in D:
                        raise ValueError("q=%s is included for g=% but its divisor %s is not" % (q, g, d))
                for gg in range(1,g):
                    if q not in gq_dict[gg]:
                        raise ValueError("g=%s is included for q=%s but g=%s is not" % (g, q, gg))
                self.gq.append((g, q))
        # We want to do low q, high g first.  This way we have all the values for a given q,
        # but do the high dimension cases first since these will take the longest.
        self.gq.sort(key=lambda pair: (pair[1], -pair[0]))

    class StageGenerateSimple(Stage):
        name = 'Generate Simple'
        shortname = 'GenSimp'
        class Task(GenericTask):
            def run(self, db=None):
                stage = self.stage
                controller = stage.controller
                g, q = self.g, self.q
                def make_simples():
                    R = ZZ['x']
                    for Ppoly in self.enumerate(WeilPolynomials(2*g, q)):
                        Lpoly = R(Ppoly.reverse())
                        IC = IsogenyClass(Lpoly=Lpoly, db=db)
                        try:
                            invs, mult = IC.simplepow_brauer_data
                        except ValueError:
                            continue
                        if mult != 1:
                            continue
                        yield IC
                self.save(make_simples())

    class StageGenerateAll(Stage):
        name = 'Generate All'
        shortname = 'GenAll'
        class Task(GenericTask):
            @lazy_attribute
            def input_data(self):
                gs = range(1, self.g+1)
                q = self.q
                return [(g, self.stage.input[0].format(g=g, q=q)) for g in gs]
            def run(self, db=None):
                stage = self.stage
                controller = stage.controller
                simples = {g: list(controller.load(filename, db=db)) for (g, filename) in self.input_data}
                self.log("Loaded")
                def make_all():
                    for split in Partitions(self.g):
                        self.log("Starting %s"%split)
                        split_mD = Counter(split)
                        it = cartesian_product_iterator([CombosWithReplacement(simples[g], c) for g, c in sorted(split_mD.items())])
                        for factors in self.enumerate(it):
                            if len(factors) == 1:
                                factors = factors[0]
                            else:
                                factors = sum(factors, ())
                            factors = Counter(factors).items()
                            factors.sort(key=lambda ICpair: (ICpair[0].g, ICpair[0].poly))
                            yield IsogenyClass.by_decomposition(factors, db=db)
                self.save(make_all())
    class StageBasechange(Stage):
        name = 'Basechange'
        shortname = 'Bchange'
        @lazy_attribute
        def tasks(self):
            return [self.Task(g, q, self) for g, q in self.controller.gq if q.is_prime()]

        class Task(GenericTask):
            def __init__(self, g, p, stage):
                self.g, self.p, self.stage = g, p, stage
                self.logheader = stage.controller.logheader.format(g=g, q=p, name=stage.shortname)
                self.qs = [q for (g, q) in stage.controller.gq if g == self.g and q%p == 0]
                self.rs = [q.is_prime_power(get_data=True)[1] for q in self.qs]

            def _logdata(self):
                return '_%s_%s' % (self.g, self.p)
            @lazy_attribute
            def input_data(self):
                return [(r, self.stage.input[0].format(g=self.g, q=q)) for q, r in zip(self.qs, self.rs)]

            @lazy_attribute
            def output_data(self):
                outputs = []
                for i, (filename, attributes) in enumerate(self.stage.output):
                    if i < 2:
                        # We use r=0 and r=-1 to encode the output files for endomorphism alg data
                        outputs.append((-i, filename.format(g=self.g, p=self.p), attributes))
                    else:
                        for r in self.rs:
                            outputs.append((r, filename.format(g=self.g, q=self.p^r), attributes))
                return outputs

            @lazy_attribute
            def donefiles(self):
                return [self._done(filename) for r, filename, attributes in self.output_data]

            def run(self, db=None):
                process = psutil.Process(os.getpid())
                stage = self.stage
                controller = stage.controller
                g, p = self.g, self.p
                self.log("Loading inputs; memory %s" % (process.memory_info()[0]))
                in_db = {r: {IC.label: IC for IC in controller.load(filename, db=db)} for (r, filename) in self.input_data}
                geometric_degrees = defaultdict(set)
                pair_lcms = defaultdict(set)
                for r, ICs in in_db.items():
                    self.log("Computing geometric extension degrees r=%s; IClen %s, memory %s" % (r, len(ICs), process.memory_info()[0]))
                    for IC in self.enumerate(ICs.values()):
                        geometric_degrees[r].add(IC.geometric_extension_degree)
                    # Need to base change to all degrees in the database,
                    # even if these aren't pair lcms
                    for s in in_db.keys():
                        if s > r and s % r == 0:
                            pair_lcms[r].add(s // r)
                    for a, b in Combinations(geometric_degrees[r], 2):
                        pair_lcms[r].add(lcm(a,b))
                self.log("pair_lcms %s" % (dict(pair_lcms)))
                base_changes = defaultdict(lambda: defaultdict(list))
                simple_factors = {}
                multiplicity_records = []
                def compute_endomorphism_algebra(extension_class, base_class, extension_degree):
                    for simple_factor, mult in extension_class.decomposition:
                        SFL = simple_factor.label
                        BCR = BaseChangeRecord(base_class.label, SFL, extension_degree, mult)
                        multiplicity_records.append(BCR)
                        if SFL not in simple_factors:
                            simple_factors[SFL] = simple_factor
                    if extension_degree == base_class.geometric_extension_degree:
                        # record some geometric data
                        base_class.geometric_center_dim = extension_class.center_dim
                        base_class.is_geometrically_simple = extension_class.is_simple
                        for m in range(1, 6):
                            setattr(base_class, 'geom_dim%s_factors', sum(mult for simple_factor, mult in extension_class.decomposition if simple_factor.g == m))
                            setattr(base_class, 'geom_dim%s_distinct', len([mult for simple_factor, mult in extension_class.decomposition if simple_factor.g == m]))
                for r, ICs in in_db.items():
                    self.log("Base changing r=%s; memory %s" % (r, process.memory_info()[0]))
                    lcms = pair_lcms[r]
                    q = p^r
                    for IC in self.enumerate(ICs.values()):
                        geom_degree = IC.geometric_extension_degree
                        for s in lcms:
                            BC = IsogenyClass(Lpoly=base_change(IC.Lpoly, s, g=g, q=q), db=db)
                            base_changes[BC][r].append(IC)
                            if geom_degree % s == 0:
                                compute_endomorphism_algebra(BC, IC, s)
                                if s == geom_degree:
                                    # record some geometric data
                                    IC.geometric_center_dim = BC.center_dim
                        for s in geom_degree.divisors():
                            if s not in lcms:
                                BC = IsogenyClass(Lpoly=base_change(IC.Lpoly, s, g=g, q=q), db=db)
                                compute_endomorphism_algebra(BC, IC, s)
                for r, ICs in in_db.items():
                    self.log("Computing primitive models r=%s; memory %s" % (r, process.memory_info()[0]))
                    for IC in self.enumerate(ICs.values()):
                        IC.twist_dict = defaultdict(dict) # filled in below
                        IC.primitive_models = None # filled in below
                        if r != 1 and IC in base_changes:
                            models = sum(base_changes[IC].values(), [])
                            models.sort(key=lambda mod: (mod.q, mod.poly))
                            IC.primitive_models = []
                            for model in models:
                                if model not in base_changes:
                                    IC.primitive_models.append(model.label)
                self.log("Computing twists; memory %s" % (process.memory_info()[0]))
                for BC, ICs in self.enumerate(base_changes.items()):
                    for r in self.rs:
                        if len(ICs[r]) > 1: # rules out r=1
                            for IC in ICs[r]:
                                for JC in ICs[r]:
                                    if IC == JC:
                                        continue
                                    new_deg = BC.r // r
                                    D = IC.twist_dict[JC.label]
                                    for old_deg in list(D):
                                        if old_deg % new_deg == 0:
                                            del D[old_deg]
                                    for old_deg in D:
                                        if new_deg % old_deg == 0:
                                            break
                                    else:
                                        D[new_deg] = BC.label
                self.log("Setting twists")
                for r, ICs in self.enumerate(in_db.items()):
                    for IC in ICs.values():
                        IC.twists = []
                        for JClabel, D in IC.twist_dict.items():
                            for deg, BClabel in D.items():
                                IC.twists.append([JClabel, BClabel, deg])
                        IC.twists.sort(key = lambda x: (x[2], map(signed_class_to_int, x[0].split('.')[2].split('_'))))
                simple_factors = sorted(simple_factors.values(), key = lambda SF: (SF.g, SF.q, SF.poly))
                for r, filename, attributes in self.output_data:
                    if r < 0:
                        self.log("Saving simple factors; %s total, memory %s" % (len(simple_factors), process.memory_info()[0]))
                        controller.save(filename, simple_factors, attributes, self.t0)
                    elif r == 0:
                        self.log("Saving multiplicity records; %s total, memory %s" % (len(multiplicity_records), process.memory_info()[0]))
                        controller.save(filename, multiplicity_records, attributes, self.t0, cls=BaseChangeRecord)
                    else:
                        self.log("Saving bc isogeny classes r=%s; %s total, memory %s" % (r, len(in_db[r]), process.memory_info()[0]))
                        controller.save(filename, in_db[r].values(), attributes, self.t0)
    class StageFixEndalgData(Stage):
        name = 'FixEndalgData'
        shortname = 'FixEndalg'
        @lazy_attribute
        def tasks(self):
            return [self.Task(g, q, self) for g, q in self.controller.gq if q.is_prime()]
        class Task(GenericTask):
            def __init__(self, g, p, stage):
                self.g, self.p, self.stage = g, p, stage
                self.kwds = {'g':g, 'p':p}
                self.logheader = stage.controller.logheader.format(g=g, q=p, name=stage.shortname)
                self.qs = [q for (g, q) in stage.controller.gq if g == self.g and q%p == 0]
                self.rs = [q.is_prime_power(get_data=True)[1] for q in self.qs]
            def _logdata(self):
                return '_%s_%s' % (self.g, self.p)
            @lazy_attribute
            def input_data(self):
                return [(None, self.stage.input[0].format(g=self.g, p=self.p))] + [(r, self.stage.input[1].format(g=self.g, q=q)) for q, r in zip(self.qs, self.rs)]
            def run(self, db=None):
                stage = self.stage
                controller = stage.controller
                g, p = self.g, self.p
                sources = [(r, controller.load(filename, db=db)) for (r, filename) in self.input_data[1:]]
                seen = set()
                def make_all():
                    for r, source in sources:
                        for IC in source:
                            if IC.label in seen:
                                continue
                            if r is not None and not IC.is_simple:
                                continue
                            seen.add(IC.label)
                            yield IC
                self.save(make_all())
    class StageCombine(Stage):
        name = 'Combine'
        shortname = 'Combo'
        class Task(GenericTask):
            def ready(self):
                # We don't demand the presence of external data files.  The corresponding columns
                # are designed to be null if not present
                return all(ope(self._done(data[-1])) for data in self.input_data[:2])
            def run(self, db=None):
                stage = self.stage
                controller = stage.controller
                sources = [controller.load(filename, db=db, allow_empty=True) for (_, filename) in self.input_data]
                data = defaultdict(list)
                for source in sources:
                    for IC in source:
                        # We want to sort by poly since g and q are fixed
                        data[tuple(IC.poly)].append(IC)
                self.log("Data loaded")
                def make_all():
                    for key in self.enumerate(sorted(data.keys())):
                        ICs = data[key]
                        yield IsogenyClass.combine(ICs, db=db)
                self.save(make_all())
    class StageCollate(Stage):
        name = 'Collate'
        shortname = 'Coll'
        @lazy_attribute
        def tasks(self):
            return [self.Task(self)]
        class Task(GenericTask):
            def __init__(self, stage):
                self.stage = stage
                self.logheader = stage.controller.logheader.format(g=0, q=0, name=stage.shortname)
            def _logdata(self):
                return ''
            @lazy_attribute
            def input_data(self):
                gqs = sorted(self.stage.controller.gq)
                gps = [(g, q) for g,q in gqs if q.is_prime()]
                pair_lists = [[{'g':g, 'q':q} for g,q in gqs], [{'g':g, 'p':p} for g,p in gps], [{'g':g, 'p':p} for g,p in gps]]
                return [[((g,q), ifile.format(**kwds)) for kwds in pair_list] for ifile, pair_list in zip(self.stage.input, pair_lists)]
            def ready(self):
                return all(ope(self._done(data[1])) for table in self.input_data for data in table)
            @lazy_attribute
            def donefiles(self):
                return [self._done(output) for output, attributes in self.stage.output]
            def run(self, db=None):
                stage = self.stage
                controller = stage.controller
                # We need to rewrite some headers (a temporary fix)
                #headers = ['text:smallint:integer:integer[]:text:text[]:smallint:double precision[]:smallint:text[]:text[]:smallint:numeric[]:text:numeric:numeric[]:text:integer:smallint:smallint:smallint:smallint:text[]:boolean:jsonb:boolean:text[]:smallint[]:text[]:smallint:smallint:smallint:smallint:smallint:smallint:smallint:smallint:smallint:smallint:integer[]:integer:boolean:boolean:numeric:numeric[]:numeric:numeric[]:numeric\n', None, None]
                headers = [None, None, None]
                #replacers = [[], [('-','')], []]
                replacers = [[], [], []]
                # We don't use load and save since we're just copying lines to the output
                for (outfile, _), infiles, H, R in zip(stage.output, self.input_data, headers, replacers):
                    header_written = False
                    with open(outfile, 'w') as Fout:
                        for (g,q), filename in infiles:
                            print filename # should go to logfile once we set up logging appropriately
                            with open(filename) as F:
                                for i, line in enumerate(F):
                                    if i > 2:
                                        for Ri, Ro in R:
                                            line = line.replace(Ri, Ro)
                                        Fout.write(line)
                                    elif not header_written:
                                        if i == 1 and H is not None:
                                            # This is a temporary fix so that we don't need to change all of the headers
                                            Fout.write(H)
                                        else:
                                            Fout.write(line)
                                            if i == 2:
                                                header_written = True
    class StageRecompute(Stage):
        name = 'Recompute'
        shortname = 'Recomp'
        def __init__(self, controller, input, output):
            Stage.__init__(self, controller, input, output)
            self.fixes = [x.strip() for x in controller.cfgp.get('StageRecompute', 'fix').split(',')]
        class Task(GenericTask):
            def run(self, db=None):
                stage = self.stage
                controller = stage.controller
                source = controller.load(self.input_data[0][1], db=db)
                def make_all():
                    for IC in source:
                        # If we delete the attributes they will be recomputed
                        for fix in stage.fixes:
                            if hasattr(IC, fix):
                                delattr(IC, fix)
                        yield IC
                self.save(make_all())
    class StageReorderDecomposition(Stage):
        name = 'ReorderDecomposition'
        shortname = 'Reorder'
        class Task(GenericTask):
            def run(self, db=None):
                stage = self.stage
                controller = stage.controller
                sources = [controller.load(filename, db=db) for _, filename in self.input_data]
                def make_all():
                    for IC in sources[0]:
                        if not IC.is_simple:
                            # sort simple_distinct, simple_multiplicities, number_fields and galois_groups
                            decomp_data = zip(IC.decomposition, IC.simple_distinct, IC.simple_multiplicities, IC.number_fields, IC.galois_groups)
                            decomp_data.sort(key=lambda data: (data[0][0].g, data[0][0].poly))
                            IC.decomposition, IC.simple_distinct, IC.simple_multiplicities, IC.number_fields, IC.galois_groups = map(list, zip(*decomp_data))
                        yield IC
                def make_bc():
                    for IC in sources[1]:
                        if not IC.is_primitive:
                            models = [IsogenyClass(label=model) for model in IC.primitive_models]
                            models.sort(key=lambda model: (model.q, model.poly))
                            IC.primitive_models = [model.label for model in models]
                        yield IC
                self.save(make_all(), make_bc())
    class StageFixTwists(Stage):
        name = 'Fix Twists'
        shortname = 'FixTwist'
        class Task(GenericTask):
            def run(self, db=None):
                stage = self.stage
                controller = stage.controller
                source = controller.load(self.input_data[0][1], db=db)
                def make_all():
                    for IC in source:
                        twists = IC.twists
                        TD = defaultdict(dict)
                        for JC, BC, n in twists:
                            for old_deg in list(TD[JC]):
                                if old_deg % n == 0:
                                    del TD[JC][old_deg]
                            for old_deg in TD[JC]:
                                if n % old_deg == 0:
                                    break
                            else:
                                TD[JC][n] = BC
                        twists = []
                        for JC, D in TD.items():
                            for deg, BC in D.items():
                                twists.append([JC, BC, deg])
                        twists.sort(key = lambda x: (x[2], map(signed_class_to_int, x[0].split('.')[2].split('_'))))
                        IC.twists = twists
                        yield IC
                self.save(make_all())
    class StageFixCols(Stage):
        name = 'Fix Columns'
        shortname = 'FixCols'
        class Task(GenericTask):
            @lazy_attribute
            def input_data(self):
                gs = range(1, self.g+1)
                q = self.q
                return [(g, self.stage.input[0].format(g=g, q=q)) for g in gs]
            def run(self, db=None):
                stage = self.stage
                controller = stage.controller
                mids = {g: list(controller.load(filename, db=db)) for (g, filename) in self.input_data}
                lookup = {IC.label: IC for ICs in mids.values() for IC in ICs}
                final = mids[self.g]
                for IC in final:
                    IC.decomposition = [(lookup[label], mult) for label, mult in IC.decomp]
                self.save(final)
    class StageDBLoad(Stage):
        name = 'DBLoad'
        shortname = 'DBLoad'
        class Task(GenericTask):
            def run(self, db=None):
                outfile, _ = self.stage.output[0]
                outfile = outfile.format(g=self.g, q=self.q)
                cols = [(attr.__name__, attr.pg_type) for attr in IsogenyClass.__dict__.values() if isinstance(attr, PGOld)]
                col_names, col_types = zip(*cols)
                selecter = SQL("SELECT {0} FROM av_fqisog WHERE g={1} AND q = {2}").format(SQL(", ").join(map(Identifier, col_names)), Literal(self.g), Literal(self.q))
                self.log('DBloaded')
                header = "%s\n%s\n\n" % (":".join(col_names), ":".join(col_types))
                db._copy_to_select(selecter, outfile, header=header, sep=":")
                self.stage.controller._finish(outfile, self.t0)
    class StageGenerateMissing(Stage):
        name = 'GenerateMissing'
        shortname = 'GenMissing'
        class Task(GenericTask):
            @lazy_attribute
            def input_data(self):
                g, q = self.g, self.q
                p, r = q.is_prime_power(get_data=True)
                if r % 2 == 1: # Only missing data when q is not a square
                    gs = range(1, g+1)
                else:
                    gs = [g]
                return [(gg, self.stage.input[0].format(g=gg, q=q)) for gg in gs]
            def run(self, db=None):
                self.log("Starting")
                stage = self.stage
                controller = stage.controller
                sources = {g: list(controller.load(filename, db=db)) for (g, filename) in self.input_data}
                lookup = {IC.label: IC for g, ICs in sources.items() for IC in ICs}
                q = self.q
                SC = IsogenyClass(poly=[1,0,-2*q,0,q**2], db=db)
                self.log("Sources loaded")
                def make_all():
                    if self.g % 2 == 0:
                        yield IsogenyClass.by_decomposition([(SC, self.g // 2)], db=db)
                    for g, olddata in sources.items():
                        if (self.g - g) % 2 != 0:
                            # Missing simple isogeny class is in g=2
                            continue
                        self.log("Starting g=%s"%g)
                        for IC in self.enumerate(olddata):
                            IC.set_from_old() # Set the new-style attributes using the ones from the database
                            if g == self.g:
                                yield IC
                            else:
                                factors = [(lookup[label], mult) for label, mult in IC.decomp] + [(SC, (self.g - g)//2)]
                                yield IsogenyClass.by_decomposition(factors, db=db)
                self.save(make_all())

    class StageBasechangeOld(Stage):
        name = 'BasechangeOld'
        shortname = 'BCO'
        @lazy_attribute
        def tasks(self):
            return [self.Task(g, q, self) for g, q in self.controller.gq if q.is_prime()]
        class Task(GenericTask):
            def __init__(self, g, p, stage):
                self.g, self.p, self.stage = g, p, stage
                self.logheader = stage.controller.logheader.format(g=g, q=p, name=stage.shortname)
                self.qs = [q for (g, q) in stage.controller.gq if g == self.g and q%p == 0]
                self.rs = [q.is_prime_power(get_data=True)[1] for q in self.qs]
            def _logdata(self):
                return '_%s_%s' % (self.g, self.p)

            @lazy_attribute
            def input_data(self):
                return [(r, self.stage.input[0].format(g=self.g, q=q)) for q, r in zip(self.qs, self.rs)]
            @lazy_attribute
            def output_data(self):
                stage = self.stage
                filename, attributes = stage.output[0]
                return [(r, filename.format(g=self.g, q=self.p^r), attributes) for r in self.rs]

            @lazy_attribute
            def donefiles(self):
                return [filename[:-4] + '.done' for (r, filename, attributes) in self.output_data]

            def run(self, db=None):
                g, p = self.g, self.p
                self.log("Starting g=%s p=%s"%(g,p))
                stage = self.stage
                controller = stage.controller
                in_db = {r: {IC.label: IC for IC in controller.load(filename, db=db)} for (r, filename) in self.input_data}
                base_changes = defaultdict(lambda: defaultdict(list))
                for ICs in in_db.values():
                    for IC in ICs.values():
                        IC.primitive_models = []
                self.log("Loading complete")
                for r, ICs in in_db.items():
                    q = p**r
                    self.log("Starting q=%s"%(q))
                    for IC in self.enumerate(ICs.values()):
                        if IC.primitive_models: # already a base change
                            continue
                        for s in self.rs:
                            # Note that this s differs from the s in BaseChange by a factor of r
                            if s > r and s%r == 0:
                                BC = IsogenyClass(Lpoly=base_change(IC.Lpoly, s//r, g=g, q=q), db=db)
                                # Get the version of BC with everything computed
                                BC = in_db[s][BC.label]
                                BC.primitive_models.append(IC.label)
                # Check to see that we agree with the old version when present
                self.log("Checking consistency....")
                for ICs in in_db.values():
                    for IC in ICs.values():
                        # Sorting?
                        if IC.prim_models is not None and any(PM not in IC.primitive_models for PM in IC.prim_models):
                            print "Prim models mismatch (%s): %s vs %s" % (IC.label, IC.prim_models, IC.primitive_models)
                            #raise RuntimeError((IC.prim_models, IC.primitive_models))
                        IC.prim_models = IC.primitive_models
                        IC.is_prim = not IC.primitive_models
                self.log("Saving....")
                for r, filename, attributes in self.output_data:
                    controller.save(filename, in_db[r].values(), attributes, self.t0)
    class StageCheck(Stage):
        name = 'Check'
        shortname = 'Check'
        class Task(GenericTask):
            @lazy_attribute
            def input_data(self):
                return [(None, self.stage.input[0].format(g=self.g, q=self.q))]
            def run(self, db=None):
                stage = self.stage
                controller = stage.controller
                source = controller.load(self.input_data[0][1], db=db)
                db_data = {rec['label']: rec for rec in db.av_fqisog.search({'g':self.g, 'q': self.q}, sort=[])}
                reports = []
                for IC in source:
                    report = IC.validation_report(db_data.get(IC.label))
                    if report:
                        reports.append(report)
                outfile, attributes = self.stage.output[0]
                outfile = outfile.format(g=self.g, q=self.q)
                if reports:
                    controller.save(outfile, reports, attributes, self.t0, cls=ValidationReport)
                else:
                    controller._finish(outfile, self.t0)
    class StagePolyOut(Stage):
        name = 'PolyOut'
        shortname = 'PolyOut'
        class Task(GenericTask):
            def run(self, db=None):
                stage = self.stage
                controller = stage.controller
                if (self.g == 1 and self.q > 499) or (self.g == 2 and self.q > 211) or (self.g == 3 and self.q > 13):
                    infile = self.input_data[0][1]
                    self.save(controller.load(infile))

    def load(self, filename, start=None, stop=None, cls=None, old=False, db=None, allow_empty=False):
        """
        Iterates over all of the isogeny classes stored in a file.
        The data contained in the file is specified by header lines: the first giving the
        column names (which are identical to PGType attributes of a PGSaver `cls`)
        the second giving postgres types (unused here), and the third blank.

        We use ':' as a separator.

        INPUT:

        - ``filename`` -- the filename to open
        - ``start`` -- if not ``None``, the line to start from, indexed so that 0 is
            the first line of data.  Note that this will be the fourth line of the
            file because of the header.
        - ``stop`` -- if not None``, the first line not to read, indexed as for ``start``
        - ``cls`` -- the PGSaver class to load data into
        - ``old`` -- if True, uses 'old_'+attr for attributes in the header if possible
        """
        if cls is None:
            cls = self.default_cls
        def fix_attr(attr):
            if old and hasattr(cls, 'old_' + attr):
                attr = 'old_' + attr
            return attr
        if ope(filename):
            with open(filename) as F:
                for i, line in enumerate(F):
                    if i == 0:
                        header = map(fix_attr, line.strip().split(':'))
                    elif i >= 3 and (start is None or i-3 >= start) and (stop is None or i-3 < start):
                        yield cls.load(line.strip(), header, db=db)
        elif not allow_empty:
            raise ValueError("%s does not exist" % (filename))

    def _finish(self, outfile, t0):
        donefile = outfile[:-4] + '.done'
        with open(donefile, 'w') as F:
            F.write(str(datetime.utcnow() - t0)+'\n')

    def save(self, filename, isogeny_classes, attributes, t0, cls=None, force=False):
        """
        INPUT:

        - ``filename`` -- a filename to write
        - ``isogeny_classes`` -- an iterable of instances to write to the file
        - ``attributes`` -- a list of attributes to save to the file
        - ``t0`` -- the time this computation was started (for recording in the .done file)
        - ``cls`` -- the class of the entries of ``isogeny_classes``
        - ``force`` -- if True, will allow overwriting an existing file
        """
        # This should be re-enabled once no longer debugging
        # if not force and ope(filename):
        #     raise ValueError("File %s already exists"%filename)
        if cls is None:
            cls = self.default_cls
        types = [getattr(cls, attr) for attr in attributes]
        header = [':'.join(attributes),
                  ':'.join(attr.pg_type for attr in types),
                  '\n']
        with open(filename, 'w') as F:
            F.write('\n'.join(header))
            for isog in isogeny_classes:
                F.write(isog.save(attributes) + '\n')
        self._finish(filename, t0)

    def combine_all(self):
        # Manual sorting to fix the critical bug; should be replaced by a stage
        old_header = None
        type_line = None
        id_cntr = 1
        with open(opj(self.basedir, "all.txt"), 'w') as Fout:
            for g, q in sorted(self.gq):
                filename = opj(self.basedir, 'complete/weil_actual_g{g}_q{q}.txt'.format(g=g, q=q))
                with open(filename) as Fin:
                    lines = []
                    for i, line in enumerate(Fin):
                        if i == 0:
                            if old_header is None:
                                old_header = line
                                pieces = old_header.strip().split(':')
                                pieces = [piece[4:] if piece.startswith('old_') else piece for piece in pieces]
                                new_header = 'id:' + ':'.join(pieces) + '\n'
                                Fout.write(new_header)
                            elif old_header != line:
                                raise ValueError
                        elif i == 1:
                            if type_line is None:
                                type_line = line
                                Fout.write('bigint:'+line+'\n')
                            elif type_line != line:
                                raise ValueError
                        elif i == 2:
                            if line != '\n':
                                raise ValueError
                        else:
                            poly = map(int, line.split(':',4)[3][1:-1].split(','))
                            lines.append((poly, line))
                    for poly, line in sorted(lines):
                        Fout.write('%s:'%id_cntr + line)
                        id_cntr += 1

    def fix_parens(self):
        with open(opj(self.basedir, "all.txt")) as Fold:
            with open(opj(self.basedir, "fixed.txt"), 'w') as Fnew:
                for line in Fold:
                    Fnew.write(line.replace("(","[").replace(")","]"))

    def fix_ratlist(self):
        quoter = re.compile(r"(-?\d+/\d+|-?\d+)")
        with open(opj(self.basedir, "all.txt")) as Fold:
            with open(opj(self.basedir, "fixed.txt"), 'w') as Fnew:
                for i, line in enumerate(Fold):
                    if i > 2:
                        pieces = line.split(':')
                        pieces[23] = quoter.sub(r'"\1"', pieces[23])
                        line = ':'.join(pieces)
                    Fnew.write(line)

class PGSaver(object):
    @classmethod
    def load(cls, s, header, db=None):
        """
        INPUT:

        - ``s`` -- a string, giving the data defined by ``header`` as colon separated values
        - ``header`` -- a list of attribute names to fill in
        """
        data = s.split(':')
        isoclass = cls(db=db)
        for attr, val in zip(header, data):
            setattr(isoclass, attr, getattr(cls, attr)._load(val))
        return isoclass

    def save(self, header):
        """
        INPUT:

        - ``header`` -- a list of attribute names to save

        OUTPUT:

        - a string, giving a colon separated list of the desired attributes
        """
        cls = self.__class__
        return ':'.join(getattr(cls, attr)._save(getattr(self, attr)) for attr in header)

class IsogenyClass(PGSaver):
    """
    An isogeny class of abelian varieties over a finite field, as constructed from a Weil polynomial
    """
    @classmethod
    def by_decomposition(cls, factors, db=None):
        """
        INPUT:

        - ``factors`` -- a list of pairs (IC, e) where IC is a simple isogeny class and e is an exponent
        """
        Lpoly = prod(IC.Lpoly^e for IC,e in factors)
        result = cls(Lpoly=Lpoly, db=db)
        result.decomposition = factors
        result.has_decomposition = True
        return result

    @classmethod
    def combine(cls, inputs, db=None):
        """
        INPUT:

        - ``factors`` -- a list of isogeny classes.  Attributes on the resulting isogeny class
        will be set from the inputs.  If there are any mismatches, a ValueError is raised.
        """
        all_attrs = {}
        for IC in inputs:
            D = IC.__dict__
            for key, val in D.items():
                if key in all_attrs:
                    if val != all_attrs[key]:
                        raise ValueError("Two different values for %s: %s and %s"%(key, val, all_attrs[key]))
                else:
                    all_attrs[key] = val
        IC = cls(db=db)
        for key, val in all_attrs.items():
            setattr(IC, key, val)
        return IC

    def __init__(self, Lpoly=None, poly=None, label=None, db=None):
        # All None is allowed since the load method writes to the fields outside the __init__ method
        if Lpoly is not None:
            if not isinstance(Lpoly, Polynomial):
                raise TypeError("Lpoly %s has type %s"%(Lpoly, type(Lpoly)))
            self.Lpoly = Lpoly
        if poly is not None:
            self.poly = poly
        if label is not None:
            self.label = label
        #if db is None:
        #    # Since we're multiprocessing we need a new database connection for each process
        #    raise RuntimeError
        #    db = PostgresDatabase()
        self.db = db
        # One of the gotchas of lazy_attributes is that hasattr triggers the computation,
        # so we have to store the existence of the decomposition in a separate variable.
        self.has_decomposition = False

    def __eq__(self, other):
        return isinstance(other, IsogenyClass) and self.Lpoly == other.Lpoly

    def __ne__(self, other):
        return not (self == other)

    def __hash__(self):
        return hash(self.Lpoly)

    @pg_integer_list
    def poly(self):
        g, q, coeffs = self.label.split('.')
        self.g = g = ZZ(g)
        self.q = q = ZZ(q)
        coeffs = map(signed_class_to_int, coeffs.split('_'))
        coeffs = [ZZ(1)] + coeffs + [q^i * c for (i, c) in enumerate(reversed(coeffs[:-1]), 1)] + [q^g]
        return coeffs

    @lazy_attribute
    def Lpoly(self):
        return ZZ['x'](self.poly)

    @pg_text_old # FIXME once no longer porting
    def label(self):
        if 'short_curve_counts' in self.__dict__:
            # We allow specifying q and short_curve_counts
            Lpoly = self._from_curve_counts() # also sets g
        else:
            Lpoly = self.Lpoly
        g, q = self.g, self.q
        return '%s.%s.%s' % (
            g, q, '_'.join(signed_cremona_letter_code(Lpoly[i]) for i in range(1,g+1)))

    @lazy_attribute
    def Ppoly(self):
        return self.Lpoly.reverse()

    @pg_smallint_old # FIXME once no longer porting
    def g(self):
        d = self.Lpoly.degree()
        if d % 2 == 1:
            raise ValueError("Must have even degree")
        return d // 2

    @pg_integer_old # FIXME once no longer porting
    def q(self):
        return self.Ppoly[0].nth_root(self.g)

    @lazy_attribute
    def p(self):
        p, r = self.q.is_prime_power(get_data=True)
        self.r = r
        return p

    @lazy_attribute
    def r(self):
        p, r = self.q.is_prime_power(get_data=True)
        self.p = p
        return r

    @lazy_attribute
    def Ppoly_factors(self):
        return self.Ppoly.factor()

    @pg_rational_mults
    def slopes(self):
        p, r, g = self.p, self.r, self.g
        np = self.Lpoly.change_ring(Qp(p)).newton_polygon()
        try:
            return [(np(i+1) - np(i)) / r for i in range(2*g)]
        except SignError:
            print self.label
            raise

    @pg_smallint_old # FIXME once no longer porting
    def p_rank(self):
        return self.slopes.count(0)

    @pg_float8_list
    def angles(self):
        """
        A list of arguments (angles) of the roots of the L-polynomial
        that are in the upper half plane,
        divided by pi, with multiplicity, coerced into RR.
        """
        if self.has_decomposition:
            return sorted(sum((IC.angles * e for IC, e in self.decomposition), []))
        Ppoly = self.Ppoly
        q = self.q
        x = Ppoly.parent().gen()
        angles = []
        # split off the real roots for numerical stability reasons
        if self.r % 2 == 0:
            p = q.isqrt()
            real_polynomials = [(x - p, [RR(0)]), (x + p, [RR(1)])]
        else:
            real_polynomials = [(x^2 - q, [RR(0), RR(1)])]
        for rpoly, real_roots in real_polynomials:
            v = Ppoly.valuation(rpoly)
            angles.extend(real_roots * v)
            Ppoly = Ppoly // rpoly^v
        roots = Ppoly.roots(CC)
        angles.extend(sum([[RR(root.arg()/pi)]*e for (root, e) in roots if root.imag() > 0],[]))
        return sorted(angles)

    @staticmethod
    def _real_poly(poly, g, q):
        R = poly.parent()
        if g == 1:
            return R([poly[1],1])
        else:
            return R([poly[g+i] +
                      sum((-q)^j * (i+2*j)/(i+j)
                          * binomial(i+j,j) * poly[g+i+2*j]
                          for j in range(1,(g+i)//2+1))
                      for i in range(g+1)])

    @lazy_attribute
    def real_Lpoly(self):
        return self.real_Ppoly.reverse()

    @lazy_attribute
    def real_Ppoly(self):
        return self._real_poly(self.Ppoly, self.g, self.q)

    @lazy_attribute
    def real_Ppoly_factors(self):
        return self.real_Ppoly.factor()

    @pg_integer_list
    def real_poly(self):
        return self.real_Lpoly.list()

    @pg_numeric_list
    def abvar_counts(self):
        """
        A list containing the number of points of an abelian variety
        in this class, defined over F_q^i, for i = 1 to max(g,10)
        """
        L, g = self.Lpoly, self.g
        x = L.parent().gen()
        return [L.resultant(x^i-1) for i in range(1,max(g,10)+1)]

    @pg_numeric_list
    def curve_counts(self):
        """
        A list containing the number of points of a possible curve
        whose Jacobian could be in this class,
        defined over F_q^i, for i = 1 to max(g,10)
        """
        L, q, g = self.Lpoly, self.q, self.g
        prec = max(g,10)
        S = PowerSeriesRing(QQ, 'x', prec+2)
        x = S.gen()
        f = S(L)/((1-x)*(1-q*x))
        return list(f.log().derivative())[:prec]

    @pg_numeric_list
    def short_curve_counts(self):
        """
        As for curve_counts, but only of length g.  The intended use
        is so that you can specify an isogeny class by providing q
        and these counts.  Note that g is determined by the length of
        the list, so you cannot provide extra counts.
        """
        return self.curve_counts[:self.g]

    def _from_curve_counts(self):
        """
        Returns a power series from the value of ``short_curve_counts``
        whose coefficient match Lpoly up to O(x^(g+1)), and sets ``g``.

        The attributes ``q`` and ``short_curve_counts`` must be set.
        """
        q = self.q
        counts = self.short_curve_counts
        self.g = g = len(counts)
        S = PowerSeriesRing(QQ, 'x', g+1)
        x = S.gen()
        return S(counts).integral().exp() * (1-x) * (1-q*x)

    @lazy_attribute
    def _newpoints(self):
        """
        The number of points of degree d (ie, we remove the points with coefficients in a subfield).

        We only go up to g, since this is all that's required for testing Jacobians.
        """
        g = self.g
        newpoints = self.curve_counts[:g] # copy
        for d in range(1, g+1):
            cnt = newpoints[d-1]
            for i in range(2, g//d + 1):
                newpoints[d*i-1] -= cnt
        return newpoints

    @pg_integer
    def curve_count(self):
        return self.curve_counts[0]

    @pg_numeric
    def abvar_count(self):
        return self.abvar_counts[0]

    @pg_text_old # FIXME once no longer porting
    def poly_str(self):
        return ' '.join(map(str, self.poly)) + ' '

    @pg_text
    def abvar_counts_str(self):
        # Trailing space so that the searching on initial segment doesn't mismatch on an initial segment
        return ' '.join(map(str, self.abvar_counts)) + ' '

    @pg_text
    def curve_counts_str(self):
        # Trailing space so that the searching on initial segment doesn't mismatch on an initial segment
        return ' '.join(map(str, self.curve_counts)) + ' '

    @lazy_attribute
    def K(self):
        """
        Should only be called for simple abvars
        """
        factors = self.Ppoly_factors
        if len(factors) != 1:
            raise ValueError("Non-simple")
        return NumberField(factors[0][0], 'a')

    @lazy_attribute
    def primes_above_p(self):
        """
        The primes in Q(F), sorted by their slope.

        Raises a ValueError if not simple
        """
        # Can't check is_simple, since that generates an infinite loop
        # if not self.is_simple:
        #     raise ValueError("Non-simple")
        p, r = self.p, self.r
        K = self.K
        a = K.gen()
        def slope_key(v):
            return a.valuation(v) / K(p^r).valuation(v)
        return sorted(K.primes_above(p), key=slope_key)

    @lazy_attribute
    def _prime_slopes(self):
        """
        The slopes in the same order as the primes above p.

        Note that this may have a different length than the isogeny class' slopes,
        since there is no repetition from residue class degree or ramification index
        """
        K = self.K
        a = K.gen()
        p, r = self.p, self.r
        return [a.valuation(v) / K(p^r).valuation(v) for v in self.primes_above_p]

    @lazy_attribute
    def simplepow_brauer_data(self):
        """
        Computes the Brauer invariants for a power of a simple isogeny class
        and checks the conditions of the Honda-Tate theorem.

        Raises a ValueError if not a power of a simple isogeny class.

        OUTPUT:

        - ``invs`` -- a list of rational numbers giving the Brauer invariants
          for the endomorphism algebra of the simple factor
          as a division algebra over its center.
        - ``mult`` -- the multiplicity of the simple factor
        """
        q, p, r = self.q, self.p, self.r
        factors = self.Ppoly_factors
        if len(factors) > 1:
            raise ValueError("Non-simple")
        factor, power = factors[0]
        invs = []
        # Need to use the ordering on the slopes corresponding to primes_above_p, rather than self.slopes
        for v, vslope in zip(self.primes_above_p, self._prime_slopes):
            inv = vslope * v.residue_class_degree() * v.ramification_index()
            invs.append(inv - inv.floor())
        e = lcm([a.denominator() for a in invs])
        # When q is not a square, the case 1-qx^2 must be handled separately.
        if not q.is_square() and factor.degree() == 2 and factor[1] == 0 and factor[0] < 0:
            e = 2
        if power % e != 0:
            raise ValueError("Honda-Tate failure")
        if power == e:
            self.is_simple = True
        else:
            self.is_simple = False
        return invs, power // e

    @lazy_attribute
    def decomposition(self):
        # We prefer to calculate from simple_distinct and simple_multiplicities
        if not self.has_decomposition:
            # We use the has_decomposition attribute to prevent infinite recursion
            self.has_decomposition = True
            return [(IsogenyClass(label=label, db=self.db), e) for label, e in zip(self.simple_distinct, self.simple_multiplicities)]
        else:
            # Lazy attributes store their result as an attribute, so the only way we can reach this branch
            # is if the other recursed back.  So we compute the decomposition from the factorization.
            factorization = []
            for factor, power in self.Ppoly_factors:
                Lfactor = factor.reverse()^power
                IC = IsogenyClass(Lpoly=Lfactor, db=self.db)
                invs, e = IC.simplepow_brauer_data
                if e != 1:
                    Lfactor = factor.reverse()^(power//e)
                    IC = IsogenyClass(Lpoly=Lfactor, db=self.db)
                IC.is_simple = True
                factorization.append((IC, e))
            factorization.sort(key=lambda pair: (pair[0].g, pair[0].poly))
            return factorization

    # As a convenience, we define multiplication and exponentiation of isogeny classes
    # using the same operation on Lpolynomials
    def __pow__(self, e):
        return IsogenyClass(Lpoly=self.Lpoly^e, db=self.db)
    def __mul__(self, other):
        if isinstance(other, basestring):
            other = IsogenyClass(label=other, db=self.db)
        elif isinstance(other, list):
            other = IsogenyClass(poly=other, db=self.db)
        elif isinstance(other, Polynomial) and other[0] == 1:
            other = IsogenyClass(Lpoly=other, db=self.db)
        elif isinstance(other, (Integer,int)) and other == 1:
            return self
        elif not isinstance(other, IsogenyClass):
            raise ValueError("Cannot multiply by %s"%other)
        return IsogenyClass(Lpoly=self.Lpoly*other.Lpoly, db=self.db)
    def __repr__(self):
        return self.label

    @pg_boolean
    def is_simple(self):
        return len(self.decomposition) == 1 and self.decomposition[0][1] == 1

    @pg_boolean
    def is_geometrically_simple(self):
        # set by base change
        g, q, s = self.g, self.q, self.geometric_extension_degree
        return IsogenyClass(Lpoly=base_change(self.Lpoly, s, g=g, q=q), db=self.db).is_simple

    @pg_text_list
    def simple_distinct(self):
        return [IC.label for IC, e in self.decomposition]

    @pg_text_list
    def simple_factors(self):
        # A list of labels of simple factors.  Duplicated factors will have "A", "B", etc appended.
        return [IC.label + cremona_letter_code(i).upper() for IC, e in self.decomposition for i in range(e)]

    @pg_smallint_list
    def simple_multiplicities(self):
        return [e for IC, e in self.decomposition]

    @pg_rational_list
    def brauer_invariants(self):
        if not self.is_simple:
            return None
            #raise ValueError("Non-simple")
        return self.simplepow_brauer_data[0]

    @pg_rational_list
    def places(self):
        """
        A list of rational numbers, giving coordinates for the second entry
        of the two-generator representation of the ideals above p.

        Only used for simple isogeny classes
        """
        if not self.is_simple:
            return None
            #raise ValueError("Non-simple")
        places = []
        p = self.p
        for v in self.primes_above_p:
            vgen = v.gens_two()[1].list()
            d = lcm(c.denominator() for c in vgen)
            valp_d, unit_d = d.val_unit(p)
            # together with p, the following will still generate the same ideal.
            # We may be able to simplify the coefficients more in the case that valp_d > 0,
            # but it's complicated.
            vgen = [((d*c)%(p^(valp_d+1)))/p^valp_d for c in vgen]
            places.append(vgen)
        return places

    @lazy_attribute
    def _nf_data(self):
        """
        List of labels for the number fields corresponding to the irreducible factors
        """
        nfs = []
        gals = []
        R = self.Ppoly.parent()
        for poly, e in self.Ppoly_factors:
            coeffs = R(pari(poly).polredbest().polredabs()).coefficients(sparse=False)
            if len(coeffs) == 3: # quadratic, where we don't need a db lookup
                sig = 0 if coeffs[0] > 0 else 2
                disc = coeffs[1]^2 - 4*coeffs[0]
                rec = {'label':'2.%s.%s.1' % (sig, abs(disc)),
                       'degree':2,
                       'galt':1}
            else:
                rec = self.db.nf_fields.lucky({'coeffs':coeffs}, projection=['label','degree','galt'], sort=[])
            if rec is None: # not in LMFDB
                nfs.append(r'\N')
                gals.append(r'\N')
            else:
                label = rec['label']
                gal = "%sT%s" % (rec['degree'], rec['galt'])
                nfs.append(label)
                gals.append(gal)
        return nfs, gals

    @pg_text
    def number_field(self):
        # Used for caching after the first stage
        if not self.is_simple:
            raise ValueError("Non-simple")
        return self.number_fields[0]

    @pg_smallint
    def center_dim(self):
        if self.has_decomposition:
            return sum(simple_factor.center_dim for simple_factor, e in self.decomposition)
        else:
            return sum(poly.degree() for poly, e in self.Ppoly_factors)

    @pg_smallint
    def geometric_center_dim(self):
        g, q, s = self.g, self.q, self.geometric_extension_degree
        return IsogenyClass(Lpoly=base_change(self.Lpoly, s, g=g, q=q), db=self.db).center_dim

    @pg_text_list
    def number_fields(self):
        return self._nf_data[0]

    @pg_text_list
    def galois_groups(self):
        return self._nf_data[1]

    @staticmethod
    def _is_significant(rel, sprec):
        m = min(map(abs, rel))
        M = max(map(abs, rel))
        if (m+1).exact_log(2) >= sprec:
            return False
        elif (M+1).exact_log(2) >= sprec:
            raise RuntimeError("Mixed significance")
        else:
            return True

    @classmethod
    def _compute_rank(cls, numbers, sprec):
        if len(numbers) == 1:
            return 1
        rel = [ZZ(c) for c in gp.lindep(numbers)]
        if cls._is_significant(rel, sprec):
            # significant relation, so we remove one of the related numbers and recurse
            for i, c in enumerate(rel):
                if c != 0:
                    numbers.pop(i)
                    return cls._compute_rank(numbers, sprec)
        else:
            # Relation not significant, so full rank
            return len(numbers)

    @pg_smallint
    def angle_rank(self):
        sprec = 25
        prec = sprec^2
        # We can't just use angles since we need higher precision (and no multiplicities)
        roots = self.Ppoly.radical().roots(ComplexField(prec), multiplicities=False)
        angles = [z.arg()/RealField(prec)(pi) for z in roots]
        angles = [angle for angle in angles if 0 < angle < 1] + [1]
        return self._compute_rank(angles, sprec) - 1

    @pg_smallint
    def geometric_extension_degree(self):
        r"""
        The smallest degree of field ext where all endomorphism are defined
        """
        Lpoly, g, q = self.Lpoly, self.g, self.q
        x = Lpoly.parent().gen()

        square = tensor_charpoly(Lpoly,Lpoly)(x/q)
        # There seems to be a difficult-to-diagnose bug in factor, which causes
        # PariError: bug in gerepile, significant pointers lost, please report
        # Scaling square so that it has integral coefficients helps but
        # doesn't completely solve the problem.
        square *= square.denominator()

        try:
            factors = square.factor()
        except PariError:
            print "Pari error in %s" % self.label
            # Try to work around the bug using squarefree decomposition
            sqfree = square.squarefree_decomposition()
            factors = []
            for fac, e in sqfree:
                for subfac, sube in fac.factor():
                    factors.append((subfac, e*sube))
        fieldext = ZZ(1)
        for factor, power in factors:
            m = Cyclotomic().order(factor)
            if m > 0:
                fieldext = fieldext.lcm(m)

        return fieldext

    @pg_text_list
    def primitive_models(self):
        # Will be set by the base change stage code
        print self.label
        raise RuntimeError

    @pg_boolean
    def is_primitive(self):
        return self.primitive_models is None

    @pg_jsonb
    def twists(self):
        # Will be set by the base change stage code
        raise RuntimeError
    @pg_integer
    def twist_count(self):
        return len(self.twists)

    # The following columns are used in av_fq_endalg_data as synonyms

    @pg_text
    def extension_label(self):
        return self.label

    @pg_text
    def center(self):
        # Used as a synonym when storing endomorphism data for simple factors of base changes
        return self.number_field

    @pg_text
    def galois_group(self):
        # Also used for caching after the first stage
        if not self.is_simple:
            raise ValueError("Non-simple")
        return self.galois_groups[0]

    @pg_smallint
    def divalg_dim(self):
        return (2*self.g // self.center_dim)^2


    # Functions for determining whether the isogeny class corresponding to a Weil polynomial
    # has a principal polarization/a Jacobian

    # The following functions are ported from Everett Howe's code IsogenyClasses.magma

    @staticmethod
    def reduced_resultant(f, g):
        """
        The *reduced resultant* of two polynomials ``f`` and ``g`` in ``Z[x]`` is the
        characteristic of the quotient ring ``Z[x] / (f,g)``.   If ``f`` and ``g`` have a common
        factor, this will be ``0``.  Otherwise, it is the smallest positive integer that
        lies in the ideal ``(f,g)``.
        """
        d, a, b = xgcd(f.change_ring(QQ), g.change_ring(QQ))
        if d.degree() > 0:
            return 0
        return lcm(c.denominator() for c in a)

    def modified_reduced_resultant(self, h1, h2):
        """
        Suppose ``h1`` and ``h2`` are real Weil `q`-polynomials (assumed to be coprime to one
        another), with associated Frobenius elements ``pi1`` and ``pi2`` in the centers of
        their associated endomorphism rings ``R1`` and ``R2``.  Let ``pi`` be the element
        ``(pi1, pi2)`` of ``R1 x R2``, and let ``pibar = q/pi``.  The *modified reducued
        resultant* of ``h1`` and ``h2`` is the smallest positive integer `n` such that ``(0, n)``
        lies in the subring ``Z[pi,pibar]`` of ``R1 x R2``.

        Usually, this `n` is simply the reduced resultant of ``h1`` and ``h2``, but if the
        product ``h1*h2`` is divisible by ``x^2 - 4*q``, we can sometimes divide the reduced
        resultant by 2.
        """
        q = self.q
        x = h1.parent().gen()
        d, a1, a2 = xgcd(h1.change_ring(QQ), h2.change_ring(QQ))
        if d.degree() > 0:
            return 0
        n = lcm(c.denominator() for c in a1)
        h = h1*h2
        H, rem = h.quo_rem(x^2-4*q)
        if rem == 0:
            splitelt = n*a1*h1 # 0 mod h1, n mod h2
            otherelt = splitelt + x*H
            g = gcd([ZZ(c) for c in otherelt])
            if g % 2 == 0:
                return n // 2
        return n

    @pg_smallint
    def has_principal_polarization(self):
        g, q, coeffs = self.g, self.q, self.poly
        if g == 1:
            return 1r
        if g == 2:
            # P-poly: x^4 = ax^3 + bx^2 + aqx + q^2
            # Howe, Maisner, Nart, Ritzenthaler
            # "Principally polarizable isogeny classes of abelian surfaces over finite fields"
            a = ZZ(coeffs[1])
            b = ZZ(coeffs[2])
            if (a^2 - b == q and b < 0 and all(p % 3 == 1 for p in b.prime_divisors())):
                return -1r
            else:
                return 1r
        elif self.is_simple:
            if g % 2 == 1:
                # Every odd-dimensional simple isogeny class has a principal polarization
                return 1r
            # We need to take radicals since the Ppoly could be a power of an irreducible
            plus_poly = self.real_Ppoly.radical()
            # Look at the CM field K and its real subfield K+.
            # If K/K+ is ramified at a finite prime, or if there is a prime of K+
            # that divides F - V and that is inert in K/K+,
            # then there's a principally polarized variety in the isogeny class.
            poly = self.Ppoly.radical()
            try:
                K.<F> = NumberField(poly)
                Kplus.<Fplus> = NumberField(plus_poly)
            except ValueError:
                print self.label
                print poly.factor()
                print plus_poly.factor()
                raise
            D = K.discriminant()
            Dplus = Kplus.discriminant()
            if D.abs() != Dplus^2:
                return 1r
            V = q / F
            # F -> V is complex conjugation
            try:
                conj = K.hom([V])
            except TypeError:
                print self.label
                print poly
                raise
            for PP, e in K.ideal(F - V).factor():
                # Being inert in K/K+ is the same as being equal to the complex conjugate
                if PP == PP.apply_morphism(conj):
                    return 1r
            if coeffs[g].gcd(q) == 1:
                # Otherwise, if ordinary, you can do the following:
                # Let N be the positive square root of Norm_{K/Q} (F - V).
                # (If we are in this case, the norm is a square.)
                # If q > 2, there is a principally-polarized variety in the isogeny class iff
                # N is congruent modulo q to the coefficient of x^g in the Weil polynomial.
                # If q = 2, there is a PPAV in the isogeny class if and only if N is congruent
                # modulo 4 to the coefficient of x^g in the Weil polynomial.
                Nsquared = ZZ((F-V).norm())
                N = Nsquared.isqrt()
                if Nsquared != N^2:
                    raise RuntimeError
                qq = q if q > 2 else 4
                return 1r if (N - coeffs[g]) % qq == 0 else -1r
        elif all(IC.has_principal_polarization for IC, mult in self.decomposition):
            # If every factor can be principally polarized, then so can the product
            # The converse isn't true
            return 1r
        return 0r

    @pg_smallint
    def has_jacobian(self):
        g, q, p, r = self.g, self.q, self.p, self.r
        p_rank = self.p_rank
        coeffs = self.poly
        if g == 1:
            return 1r
        elif g == 2:
            # Howe, Nart, Ritzenthaler
            # "Jacobians in isogeny classes of abelian surfaces over finite fields"
            if self.is_simple:
                # P-poly: x^4 + ax^3 + bx^2 + aqx + q^2
                a = ZZ(coeffs[1])
                b = ZZ(coeffs[2])
                if (a^2 - b == q and b < 0 and all(p % 3 == 1 for p in b.prime_divisors())
                    or p_rank == 2 and a == 0 and (b == 1-2*q or p > 2 and b == 2-2*q)
                    or p_rank == 0 and a == 0 and b == -q and (p % 12 == 11 and r % 2 == 0 or
                                                               p == 3 and r % 2 == 0 or
                                                               p == 2 and r % 2 == 1)
                    or p_rank == 0 and a == 0 and b == -2*q and (q == 2 or q == 3)):
                    return -1r
                else:
                    return 1r
            else:
                # P-poly: (x^2 - sx + q)(x^2 - tx + q) with |s| >= |t|
                if len(self.decomposition) == 1:
                    # square of an elliptic curve
                    SIC, _ = self.decomposition[0]
                    s = t = -SIC.poly[1]
                else:
                    (SIC, _), (TIC, _) = self.decomposition
                    s = -ZZ(SIC.poly[1])
                    t = -ZZ(TIC.poly[1])
                    if abs(t) > abs(s):
                        s, t = t, s
                if (abs(s - t) == 1
                    or p_rank == 2 and (s == t and t^2 - 4*q in [-3, -4, -7] or
                                        q == 2 and abs(s) == abs(t) == 1 and s != t)
                    or p_rank == 1 and r % 2 == 0 and s^2 == 4*q and (s-t).is_squarefree()
                    or p_rank == 0 and (p > 3 and abs(s) != abs(t) or
                                        p == 3 and r % 2 == 1 and s^2 == t^2 == 3*q or
                                        p == 3 and r % 2 == 0 and (s - t) % (3*p^(r//2)) != 0 or
                                        p == 2 and (s^2 - t^2) % (2*q) != 0 or
                                        q in [2,3] and s == t or
                                        q in [4,9] and s^2 == t^2 == 4*q)):
                    return -1r
                else:
                    return 1r
        elif (self.has_principal_polarization == -1 or
              self._nojac_pointcounts() or
              self._nojac_stohr_voloch() or
              self._nojac_beauville_zaytsev() or
              self._nojac_korchmaros_torres() or
              self._nojac_howe_lauter() or
              self._nojac_serre()):
            return -1r
        return 0r

    def _nojac_pointcounts(self):
        """
        Returns True if we can rule out the presence of a Jacobian by examining the point counts
        of the virtual curve.

        Namely, if the virtual curve has a negative number of degree d points for any d up to g,
        or if the number of points drops from F_{q^m} to F_{q^n} when m divides n,
        this isogeny class cannot contain a Jacobian.
        """
        counts = [None] + self.curve_counts # counts[m] = count over F_{q^m}
        for m in range(1, len(counts)):
            cnt = counts[m]
            if cnt < 0:
                return True
            for n in range(2*m, len(counts), m):
                ext_cnt = counts[n]
                if ext_cnt < cnt:
                    return True
        return False

    # The following _nojac methods are from Howe's http://ewhowe.com/Magma/IsogenyClasses.magma

    def _nojac_stohr_voloch(self):
        # For certain values of q and g, the paper [Stohr and Voloch 1986]
        # gives bounds on N_q(g) that may be better than the Oesterle bound.
        g, q, p = self.g, self.q, self.p
        if g < 3 or p < 2*g - 3:
            # No info in this case
            return
        N = self.curve_count

        # See [Stohr and Voloch 1986], Proposition 3.2, p. 15.  The proposition
        # gives a bound for non-hyperellipic curves, but the bound holds for
        # hyperelliptic curves as well, because it is greater than 2*q + 2.
        SVbound1 = g*(g-2) + (q*(2*g-3)) // (g-2)
        if N > SVbound1:
            return True
        if p < 2*g - 1:
            return

        # See [Stohr and Voloch 1986], Proposition 3.1, p. 15, together with the
        # comment following the proposition.
        SVbound2 = g*(g-1) + 2*q
        if N > SVbound2:
            return True

    def _nojac_beauville_zaytsev(self):
        # There are restrictions on when a curve can have Jacobian isogenous
        # to E^g, for certain elliptic curves E / F_q and genera g.  This
        # procedure looks at the cases where E has trace t and where t^2 - 4*q
        # lies in {-3, -4, -7, -8, -11, -19}.

        # The case of discriminants -3 and -4 is covered by an argument
        # attributed to Beauville; see [Serre 1985], pp. Se13--Se15.
        # The result is that no curve of genus greater than 1 has Weil
        # polynomial (x^2 - t*x + q)^g if t^2 - 4*q is -3 or -4.

        # The case of discriminant -7 can only happen when q is a power of 2.
        # The result is that any curve with Jacobian isogenous to E^g must be
        # the base extension of a curve over GF(2).  If there are any issues
        # with this that prevent such a curve from existing, they will be found
        # by the routine is_eliminated_by_Galois_descent.  So here we will do
        # nothing.

        # The case of discriminant -8 is considered in [Zaytsev 2014], but there
        # are mistakes in the proofs.  Zaytsev's result is that there is no
        # curve with Weil polynomial (x^2 - t*x + q)^g for 2 < g < 8 if
        # t^2 - 4*q = -8.   Note that this means that t = 2*s and q = s^2 + 2
        # for some odd integer s.  The prime divisor p of q is therefore either
        # equal to 3 or greater than 7.  The case q = 3^n can be analyzed using
        # Galois descent, so we assume the characteristic is at least 11.
        # The elliptic curve of trace t has j = 8000 and has  q + 1 - t
        # = s^2 - 2*s + 3 rational points.  Since s is odd, this quantity is
        # congruent to 2 modulo 4.  Since the characteristic is at least 11,
        # 8000 is different from 0 and from 1728 in GF(q), so Aut E = {1,-1}.
        # Then the automorphism group of the genus-1 curve E is isomorphic to
        # E(F_q) semidirect {1,-1}, and the 2-part of this is (Z/2) x (Z/2).

        # Zaytsev's argument for genus 3 is correct.
        # His argument for genus 4 is incorrect, but can be fixed.  The error
        # is in assuming that the genus-1 curve E has automorphism group {1,-1}.
        # But the 2-part of the automorphism group has order 4, so the argument
        # works.
        # The arguments for genus 5 and genus 6 work.
        # For genus 7, we argue (as does Zaytsev in some cases) that the
        # putative curve has automorphism group with Sylow 2-group large,
        # and that there is a non-central involution coming from a degree-2
        # map to E.  (Non-central because E has small 2-part of automorphism
        # group.)  Throwing in a central involution, we get a V4 subgroup.
        # What can the three quotients be?  There's E, and then two other
        # curves, of genus at most 2.  This contradicts the Kani-Rosen
        # decomposition result.

        # The case of discriminant -11 is considered in [Zaytsev 2016], but
        # again there is an error caused by missing the distinction between
        # the automorphism group of E as a curve versus as an elliptic curve.

        # For genus 3, no curve with Jacobian isomorphic to E^3 exists,
        # by an argument on Hermitian forms.

        # For genus 4 we obtain...
        # Theorem:
        # Suppose E is an elliptic curve over F_q with trace t such that
        # t^2 - 4*q = -11.  If the prime divisor of q is neither 3 nor 5 then
        # there is no genus-4 curve over F_q with Jacobian isogenous to E^4.

        # The proof involves analysing the automorphism groups of Hermitian forms,
        # as Zaytsev does.  What we can show is that a genus-4 curve with
        # Jacobian E^4 must have two commuting involutions that fit into a
        # V4 diagram whose middle quotients are E, E, and the unique genus-2
        # curve with Jacobian isomorphic to E^2.  By considering an explicit
        # model of this genus-2 curve, we find that such a diagram can only
        # exist in characteristics 3 and 5.

        # For genus 5, we again find V4 subgroups of the automorphism group,
        # with intermediate quotients of genus 1, 1, and 3.  But no such curve
        # of genus 3 exists.

        # The case of discriminant -19 is very finicky, and I have not verified
        # results, so I will not include this case here yet.

        g, q, p = self.g, self.q, self.p
        t = -self.real_Ppoly[g-1] // g
        x = polygen(ZZ)
        if self.real_Ppoly != (x-t)^g:
            return
        disc = t^2 - 4*q

        if disc in [-3, -4]:
            return (g > 1)
        elif disc == -8:
            return (3 <= g <= 7)
        elif disc == -11:
            return ((g in [3,5]) or (p not in [3,5] and g == 4))

    def _nojac_korchmaros_torres(self):
        # If an isogeny class is maximal, we can check to see whether the
        # conditions of [Korchmaros and Torres 2002] are satisfied.
        g, p, r = self.g, self.p, self.r
        if r % 2 != 0:
            return
        Q = p^(r//2)
        m = 2*Q
        x = polygen(ZZ)
        if self.real_Ppoly != (x + m)^g:
            return

        # See Corollary 1, p. 595, of [Korchmaros and Torres 2002]
        if (g > ((Q^2 - Q + 4) // 6) and
            g != (Q - 1)^2 // 4 and
            g != Q*(Q-1) // 2):
            return True

    def _nojac_howe_lauter(self):
        # Suppose q is a square, say q = s^2 with s positive or negative, 
        # and we can write the real Weil polynomial as (x - 2*s)^n * h0 for 
        # some ordinary h0 (!= 1).  If h0(2*s) is squarefree, then there is no
        # nontrivial self-dual group scheme that can be embedded in both
        # a variety with real Weil polynomial h0 and a variety with real
        # Weil polynomial (x - 2*s)^n.

        # This generalizes Corollary 12 of [Howe and Lauter 2003].
        g, q = self.g, self.q
        x = polygen(ZZ)
        h = self.real_Ppoly
        ss_factor = h.gcd(x^2 - 4*q)
        if ss_factor.degree() != 1:
            return
        s = -ss_factor[0] // 2
        # So the supersingular factor is a power of (x - 2*s)

        h0 = h
        while h0(2*s) == 0:
            h0 = h0 // (x - 2*s)
        # The whole variety can't consist of the supersingular part...
        if h0.degree() == 0:
            return

        critical_value = h0(2*s)
        return (q.gcd(critical_value) == 1 and critical_value.is_squarefree())

    def _nojac_serre(self):
        # Check to see whether an isogeny class is eliminated by the
        # `resultant = 1` argument.

        # Compute the factorization of the real Weil polynomial, the matrix of
        # pairwise modified reduced resultants, and the list of modified
        # reduced resultants for all possible splittings of the real Weil
        # polynomial.

        # Here is how we enumerate the possible splittings, if there are n prime
        # divisors:

        # For every integer i from 0 to 2^(n-1)-2, consider the binary
        # representation of 2^(n-1) + i.  This will be a string containing
        # exactly n bits, with the first one equal to 1, and with at least one 0.
        # This sequence determines the choice of factors.

        # Note: The integer corresponding to splitting the first factor off from
        # the rest is 0.  The integer corresponding to splitting the i-th factor
        # off from the rest, with i > 1, is 2^(n-1) - 1 - 2^(n-i).

        irred_factors = [f for f, e in self.real_Ppoly_factors]
        n = len(irred_factors)
        for i in srange(2^(n-1)-1):
            h0 = h1 = 1
            for b, factor in zip((2^(n-1)+i).bits(), irred_factors):
                if b == 0:
                    h0 *= factor
                else:
                    h1 *= factor
            res = self.modified_reduced_resultant(h0, h1)
            if res == 1:
                return True

    # The following columns are computed by Stefano's code
    @pg_boolean
    def zfv_is_bass(self):
        return None
    @pg_boolean
    def zfv_is_maximal(self):
        return None
    @pg_numeric
    def zfv_index(self):
        return None
    @pg_numeric_list
    def zfv_index_factorization(self):
        return None
    @pg_numeric
    def zfv_plus_index(self):
        return None
    @pg_numeric_list
    def zfv_plus_index_factorization(self):
        return None
    @pg_numeric
    def zfv_plus_norm(self):
        return None
    @pg_integer
    def size(self):
        return None
    @pg_integer
    def ppav_count(self):
        return None
    # The following columns are computed by Drew's code
    @pg_integer
    def jacobian_count(self):
        if self.g == 2 and self.q <= 211:
            # If not set, then we didn't find any Jacobians and the value is zero
            return 0
        else:
            return None
    @pg_integer
    def hyp_count(self):
        if self.g <= 2:
            # all curves are hyperelliptic
            return self.jacobian_count
        elif self.g == 3 and self.q in [2,3,5,7,9,11,13]:
            # If not set, then we didn't find any Jacobians and the value is zero
            return 0
        else:
            return None
    @pg_text_list
    def curves(self):
        if (self.g == 2 and self.q < 100) or (self.g == 3 and self.q < 7):
            # If not set, no Jacobians
            return []
        else:
            return None

    # The following functions are not used in the enumerations from the WeilPolynomial iterator,
    # but may be useful if you want to create isogeny classes by hand

    def check_weil_poly(self, prec=30):
        """
        Checks that the polynomial is a valid Weil polynomial.

        This function is not called manually since isogeny classes
        are usually constructed from the WeilPolynomial iterator,
        which makes this check redundant.

        Raises a ``ValueError`` on failure.
        """
        g, q, Lpoly, Ppoly = self.g, self.q, self.Lpoly, self.Ppoly
        if Lpoly[0] != 1:
            raise ValueError("Must have constant coefficient 1")
        for i in range(g):
            if Lpoly[2*g-i] != Lpoly[i]*q^(g-i):
                raise ValueError("Symmetry condition failed")
        if not all((z.abs()^2 - q).abs() < 2^-prec for z in Ppoly.roots(CC)):
            raise ValueError("Not all roots on circle of radius sqrt(q)")

    def check_honda_tate(self):
        """
        Checks that the Weil polynomial satisfies the conditions of the
        Honda-Tate theorem.

        This function is not called manually since isogeny classes
        are usually constructed from the WeilPolynomial iterator,
        which makes this check redundant.

        Raises a ``ValueError`` on failure.
        """
        # The following will raise a ValueError if Honda-Tate isn't satsified
        D = self.decomposition

    def validation_report(self, db_data):
        """
        Compares the data in this isogeny class with a record from the av_fqisog table of the LMFDB.

        Returns a ValidationReport object with the differences.  If no differences, then the
        returned object will evaluate to False in a boolean context.
        """
        diff_dict = {}
        for dbname, attr in self.__class__.__dict__.items():
            if dbname in ['C_cnts_str', 'pt_cnt']:
                # checked by C_cnts
                continue
            if not isinstance(attr, PGOld):
                continue
            if dbname.startswith('old_'):
                dbname = dbname[4:]
            icname, switch = attr.new_name, attr.switch
            if icname and switch: # Some old columns don't have direct new analogues
                dbval = switch(db_data.get(dbname), self.db)
                icval = getattr(self, icname)
                if icname == 'slopes':
                    dbval = [QQ(str(x).translate(None, letters)) for x in dbval]
                if dbname in ['simple_factors', 'simple_distinct']:
                    dbval = sorted(dbval)
                    icval = sorted(icval)
                if icname == 'angles':
                    max_diff = 0
                    for a in icval:
                        min_diff = 1
                        for b in dbval:
                            if abs(a-b) < min_diff:
                                min_diff = abs(a-b)
                        if min_diff > max_diff:
                            max_diff = min_diff
                    if max_diff > 1e-14:
                        diff_dict[dbname, icname] = (dbval, icval)
                elif dbname == 'poly_str':
                    if not icval in [dbval, dbval + ' ']:
                        diff_dict[dbname, icname] = (dbval, icval)
                elif dbname == 'C_cnts':
                    if icval != dbval:
                        moddedic = [x for x in icval if x != 0]
                        moddeddb = dbval[:len(moddedic)]
                        if moddedic != moddeddb:
                            diff_dict[dbname, icname] = (dbval, icval)
                elif dbname in ['is_pp', 'is_jac']:
                    if dbval != 0 and dbval != icval:
                        # See if it's been fixed
                        #if not (dbval == -1 and icval == 0 and self._nojac_pointcounts()):
                        diff_dict[dbname, icname] = (dbval, icval)
                elif dbval != icval:
                    diff_dict[dbname, icname] = (dbval, icval)
        return ValidationReport(self.label, diff_dict)

    # The following methods and attributes support access to the old schema

    def set_from_old(self, check=False):
        # Have to set the values first so that computations work
        if check:
            self.set_from_old(check=False)
        for aname, attr in self.__class__.__dict__.items():
            if not isinstance(attr, PGOld):
                continue
            newname, switch = attr.new_name, attr.switch
            if newname and switch: # Some old columns don't have direct new analogues
                oldval = getattr(self, aname)
                newval = switch(oldval, self.db)
                if check:
                    compval = getattr(self, newname)
                    if newval != compval:
                        print "%s (%s -> %s): %s vs %s"%(self.label, aname, newname, newval, compval)
                else:
                    setattr(self, newname, newval)

    @pg_jsonb_old
    def old_poly(self):
        return self.poly

    @pg_jsonb_old
    def old_angles(self):
        return self.angles

    @pg_smallint_old(new_name='angle_rank')
    def ang_rank(self):
        return self.angle_rank

    @pg_rational_mults_old(new_name='slopes')
    def slps(self):
        return self.slopes

    @pg_jsonb_old(new_name='abvar_counts')
    def A_cnts(self):
        return self.abvar_counts

    @pg_text_old(new_name='abvar_counts_str')
    def A_cnts_str(self):
        return self.abvar_counts_str

    @pg_jsonb_old(new_name='curve_counts')
    def C_cnts(self):
        return self.curve_counts

    @pg_text_old(new_name='curve_counts_str')
    def C_cnts_str(self):
        return self.curve_counts_str

    @pg_integer_old(new_name='curve_count')
    def pt_cnt(self):
        return self.curve_count

    @pg_smallint_old(new_name='has_jacobian')
    def is_jac(self):
        return self.has_jacobian

    @pg_smallint_old(new_name='has_principal_polarization')
    def is_pp(self):
        return self.has_principal_polarization

    @pg_jsonb_old(new_name='decomposition', switch=lambda x,db:[(IsogenyClass(label=label, db=db), e) for label, e in x])
    def decomp(self):
        return [[IC.label, e] for IC, e in self.decomposition]

    @pg_boolean_old(new_name='is_simple')
    def is_simp(self):
        return self.is_simple

    @pg_jsonb_old
    def old_simple_factors(self):
        return self.simple_factors

    @pg_jsonb_old
    def old_simple_distinct(self):
        return self.simple_distinct

    @pg_text_old(new_name='brauer_invariants', switch=False)#lambda x: x.split())
    def brauer_invs(self):
        return " ".join(" ".join(map(str, IC.brauer_invariants)) for IC, mult in self.decomposition)

    @pg_rational_list_old(switch=False)
    def old_places(self):
        return [IC.places for IC, mult in self.decomposition]

    @pg_jsonb_old(new_name='primitive_models', switch=False)
    def prim_models(self):
        # These are set in BasechangeOld.  Since we need to store these in the previous
        # stage, we return None by default.
        return None

    @pg_boolean_old(new_name='is_primitive', switch=False)
    def is_prim(self):
        # These are set in BasechangeOld.  Since we need to store these in the previous
        # stage, we return None by default.
        return None

    @pg_text_old(new_name='number_fields', switch=False)
    def nf(self):
        if len(self.decomposition) > 1:
            return None
        else:
            return self.number_fields[0]

    @pg_integer_old(new_name='galois_groups', switch=False)
    def galois_t(self):
        if len(self.decomposition) > 1:
            return None
        else:
            return int(self.galois_groups[0].split('T')[1])

    @pg_smallint_old(new_name=False)
    def galois_n(self):
        if len(self.decomposition) > 1:
            return None
        else:
            return int(self.galois_groups[0].split('T')[0])

for d in range(1,6):
    def dim_factors(self, d=d): # Need to bind d before it changes
        # Have to make a local copy of d since d is changing in the loop
        return sum([e for (simp_label, e) in zip(self.simple_distinct, self.simple_multiplicities)
                    if simp_label.split('.')[0] == str(d)])
    def dim_distinct(self, d=d):
        return len([simp_label for simp_label in self.simple_distinct if simp_label.split('.')[0] == str(d)])
    def geom_dim_factors(self, d=d):
        raise RuntimeError
    def geom_dim_distinct(self, d=d):
        raise RuntimeError
    # FIXME once no longer porting
    dim_factors.__name__ = 'dim%d_factors'%d
    dim_distinct.__name__ = 'dim%d_distinct'%d
    geom_dim_factors.__name__ = 'geom_dim%d_factors'%d
    geom_dim_distinct.__name__ = 'geom_dim%d_distinct'%d
    setattr(IsogenyClass, 'dim%d_factors'%d, pg_smallint_old(dim_factors))
    setattr(IsogenyClass, 'dim%d_distinct'%d, pg_smallint_old(dim_distinct))
    setattr(IsogenyClass, 'geom_dim%d_factors'%d, pg_smallint_old(geom_dim_factors))
    setattr(IsogenyClass, 'geom_dim%d_distinct'%d, pg_smallint_old(geom_dim_distinct))

class BaseChangeRecord(PGSaver):
    """
    This class is used to store rows for the `av_fq_endalg_factors` table.

    All attributes are set in the init method: the pg_* decorators are present
    just to explain how to save data to disk.
    """
    def __init__(self, base_label, extension_label, extension_degree, multiplicity):
        self.base_label = base_label
        self.extension_label = extension_label
        self.extension_degree = extension_degree
        self.multiplicity = multiplicity
    def _dummy(self):
        raise RuntimeError
    base_label = pg_text(_dummy)
    extension_label = pg_text(_dummy)
    extension_degree = pg_smallint(_dummy)
    multiplicity = pg_smallint(_dummy)

class ValidationReport(PGSaver):
    """
    This class is used to record discrepancies between a computed isogeny class
    and an isogeny class from the database.
    """
    def __init__(self, label, diff_dict):
        self.label = label
        self.diff_dict = diff_dict
    def __nonzero__(self):
        return bool(self.diff_dict)
    __bool__ = __nonzero__
    @pg_text
    def label(self):
        raise RuntimeError
    @pg_text
    def report(self):
        errors = []
        for (db_key, ic_key), (db_val, ic_val) in self.diff_dict.items():
            errors.append("%s (DB %s) <> %s (IC %s)" % (db_val, db_key, ic_val, ic_key))
        return "    ".join(errors)
    def _dummy(self):
        raise RuntimeError
    label = pg_text(_dummy)

def combine_polys():
    basedir = '/scratch/importing/avfq/polys/'
    for g in [1,2,3]:
        with open(opj(basedir, 'weil_poly_g%s.txt'%g), 'w') as Fout:
            for q in [16,25,243,256,343,512,625,729,1024]:
                infile = opj(basedir, 'weil_poly_g%s_q%s.txt' % (g, q))
                if ope(infile):
                    with open(infile) as Fin:
                        for i, line in enumerate(Fin):
                            if i > 2:
                                Fout.write(line)

def adjust_isom_data():
    basedir = '/home/roed/avfq/external/isom_data'
    cols = 'label:zfv_is_bass:zfv_is_maximal:zfv_index:zfv_index_factorization:zfv_plus_index:zfv_plus_index_factorization:zfv_plus_norm'
    typs = 'text:boolean:boolean:numeric:numeric[]:numeric:numeric[]:numeric'
    n = len(cols)
    for filename in os.listdir(basedir):
        g, q = map(ZZ, filename.split('_')[:2])
        outfilename = opj(basedir, 'weil_isomdata_g{g}_q{q}.txt'.format(g=g, q=q))
        with open(opj(basedir, filename)) as Fin:
            with open(opj(basedir, outfilename), 'w') as Fout:
                for i, line in enumerate(Fin):
                    if i == 0:
                        # There was a missing newline in the initial data, and we also need to add types
                        #if not line.startswith(cols) and line[n].isdigit():
                        #    raise ValueError("Update adjusting code")
                        Fout.write(cols+'\n')
                        Fout.write(typs+'\n\n')
                        Fout.write(line[n:])
                    else:
                        Fout.write(line)

def adjust_jac_data():
    indir = '/home/roed/'
    outdir = '/home/roed/avfq/external/jac_data/'
    cols = 'g:q:short_curve_counts:jacobian_count'
    typs = 'smallint:integer:numeric[]:integer'
    curgq = None
    curfile = None
    try:
        for g in [1,2]:
            with open(opj(indir, 'g%sjaccnts.txt'%g)) as Fin:
                for line in Fin:
                    pieces = line.split(':')
                    q = ZZ(pieces[0])
                    cnt = ZZ(pieces[-1])
                    # replace [] with {}
                    pieces[1] = pieces[1][1:]
                    pieces[-2] = pieces[-2][:-1]
                    curve_counts = '{' + ','.join(pieces[1:-1]) + '}'
                    if curgq != (g,q):
                        if curfile is not None:
                            curfile.close()
                        curgq = (g,q)
                        curfile = open(opj(outdir, 'weil_jacdata_g%s_q%s.txt'%(g,q)),'a')
                        if curfile.tell() == 0: # new file, normal case
                            curfile.write('%s\n%s\n\n' % (cols, typs))
                    curfile.write('%s:%s:%s:%s\n' % (g, q, curve_counts, cnt))
    finally:
        if curfile is not None and not curfile.closed:
            curfile.close()

def extract_isog_sizes():
    outfile = '/home/roed/avfq/external/isom_data/isog_sizes_g{g}_q{q}.txt'
    outs = {}
    max_q = {1:499, 2:128, 3:9, 4:4, 5:2, 6:0}
    try:
        for infile, include_max in [
                ('/home/stmar/isogeny_classes_easy_new.txt', True),
                ('/home/stmar/isogeny_classes_computations_ord.txt', False)]:
            with open(infile) as Fin:
                for i, line in enumerate(Fin):
                    if i and i%10000 == 0:
                        print infile, i
                    label, order_is_bass, order_is_maximal, size = line.strip().split(':')
                    if not include_max and order_is_maximal == '1':
                        # Don't add duplicates
                        continue
                    g, q, poly = label.split('.')
                    g, q = int(g), int(q)
                    if q <= max_q[g]:
                        if (g, q) not in outs:
                            outs[g,q] = open(outfile.format(g=g, q=q), 'w')
                            outs[g,q].write('label:size\ntext:integer\n\n')
                        outs[g,q].write('%s:%s\n' % (label, size))
    finally:
        print sorted(outs.keys())
        for ofile in outs.values():
            if not ofile.closed:
                ofile.close()

def split_easy_isoms():
    outfile = '/home/roed/isom_data/isomorphism_classes_easy_new_g{g}_q{q}.txt'
    outs = {}
    bad_context = []
    prev_line = None
    bad_entry = None
    try:
        with open('/home/stmar/isomorphism_classes_easy_new.txt') as Fin:
            for i, line in enumerate(Fin):
                if i and i%10000 == 0:
                    print i
                if bad_entry is not None:
                    bad_context.append((i, bad_entry[0], bad_entry[1], line.strip()))
                    bad_entry = None
                pieces = line.split(':')
                if len(pieces) != 5:
                    bad_entry = (prev_line.strip(), line.strip())
                    continue
                label = pieces[0]
                gq, poly = label.rsplit('.',1)
                ofile = outs.get(gq)
                if ofile is None:
                    g,q = gq.split('.')
                    ofile = outs[gq] = open(outfile.format(g=g, q=q), 'w')
                ofile.write(line)
                prev_line = line
        return bad_context
    finally:
        for ofile in outs.values():
            if not ofile.closed:
                ofile.close()

def show_ratios():
    from subprocess import check_output
    max_q = {1:499, 2:128, 3:9, 4:4, 5:2}
    isoms = defaultdict(int)
    totals = {}
    for g in range(1,6):
        for q in srange(2, 500):
            if q <= max_q[g] and q.is_prime_power():
                filename = '/home/roed/avfq/external/isom_data/isog_sizes_g{g}_q{q}.txt'.format(g=g,q=q)
                if ope(filename):
                    isoms[g,q] = int(check_output(['wc',filename]).split()[0])-3
                totals[g,q] = db.av_fqisog.count({'g':g, 'q':q})
    for g in range(1,6):
        for q in srange(2, 500):
            if q <= max_q[g] and q.is_prime_power():
                print '%6s/%-6s' % (isoms[g,q], totals[g,q]),
        print ''

def fix_138():
    inbase = '/home/stmar/138_output'
    outbase = '/home/roed/142_output_isog'
    for filename in os.listdir(inbase):
        with open(opj(inbase, filename)) as Fin:
            with open(opj(outbase, filename), 'w') as Fout:
                for i, line in enumerate(Fin):
                    if i > 0:
                        break
                    Fout.write(line.strip() + ':size\n')
def extract_bad():
    outfile = '/home/roed/isomorphism_classes_easy_new_bad.txt'
    blfile = '/home/stmar/isomorphism_classes_easy_new_bad_lines.txt'
    with open(blfile) as Fbl:
        bad_lines = set(map(lambda x: int(x)-1, Fbl.read().strip().split('\n')))
    with open('/home/stmar/isomorphism_classes_easy_new.txt') as Fin:
        with open(outfile, 'w') as Fout:
            for i, line in enumerate(Fin):
                if i in bad_lines:
                    print i
                    Fout.write(line)
