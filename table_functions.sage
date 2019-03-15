"""
Functions for generating the LMFDB data on isogeny classes of abelian varieties over FF_q

AUTHORS:
  -- (2016-05-11) Taylor Dupuy, Kiran S. Kedlaya, David Roe, Christelle Vincent
  -- (2019-02-02) David Roe, Everett Howe, added more tests for principal polarizations and Jacobians, refactored


Fields we want to populate with an example

label: "2.9.ab_d"
polynomial: ["1","-1","3","-9","81"]
angle_numbers (doubles): [0.23756..., 0.69210...]
number_field: "4.0.213413.1"
p-rank: 1
slopes: ["0","1/2","1/2","1"]
A_counts: ["75", "7125"]
C_counts: ["9", "87"]
known_jacobian (0,1,-1): 1
decomposition: ["9.2.-1._3"]
pricipally_polarizable (0,1,-1): 1
Brauer Invariants: inv_v( End(A_{FFbar_q})_{QQ} )=(v(\pi)/v(q))*[QQ(pi)_{v}: QQ(pi): v\vert p place of QQ(\pi)], these are stored as elements of QQ.
Primitive models: 
"""

######################################################################################################

from sage.databases.cremona import cremona_letter_code, class_to_int # for the make_label function
from sage.misc.lazy_attribute import lazy_attribute
import json, os, re, sys
opj, ope = os.path.join, os.path.exists
from collections import defaultdict
from itertools import combinations_with_replacement, imap, izip_longest
from ConfigParser import ConfigParser
from lmfdb import db

try:
    # Add the location of weil_polynomials.pyx to the load path
    sys.path.append(os.path.dirname(os.path.realpath(__file__)))
except NameError:
    pass
load("weil_polynomials.pyx")

######################################################################################################

# Some combinatorial utilities for enumerating non-simple isogeny classes from simple ones

def multiplicity_dict(L):
    L = list(L)
    return {a: L.count(a) for a in set(L)}

class CombosWithReplacement(object):
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
    Sym = SymmetricFunctions(QQ)
    p = Sym.powersum()
    if i == 0:
        return p.one()
    e = Sym.elementary()
    return e(p(e[i]).map_support(lambda A: Partition([r*c for c in list(A)])))

@cached_function
def basechange_transform(g, r, q):
    """
    Returns a transform that takes in a `q`-Weil L-polynomials of degree `2g` (constant coefficient 1)
    and returns a `q^r`-Weil L-polynomial of degree `2g` whose roots are the `r`th powers of the roots
    of the input.
    """
    f = [symfunc(i, r) for i in range(g+1)]
    coeffs = [b.coefficients() for b in f]
    exps = [[{a: list(elem).count(a) for a in set(elem)} for elem in sorted(b.support()) if list(elem) and max(elem) <= 2*g] for b in f]
    def bc(Lpoly):
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

def tensor_product(f, g):
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

def base_change(Lpoly, r, algorithm='sym', g = None, q = None, prec=53):
    if g is None:
        g = Lpoly.degree()
        assert g % 2 == 0
        g = g // 2
    if q is None:
        q = Lpoly.leading_coefficient().nth_root(g)
    if algorithm == 'approx':
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
            if i % r != 0 and ci != 0:
                raise RuntimeError("i = %s, ci = %s" % (i, ci))
            else:
                newf[i//r] = ZZ(ci)
        return RT(newf)
    else:
        return basechange_transform(g, r, q)(Lpoly)

class Cyclotomic(UniqueRepresentation):
    # Convention: phi(m) = n
    def __init__(self):
        self.polynomials = {}
        self.mbound = 0
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
        return 0r
    elif code.startswith('a'):
        return -class_to_int(code[1:])
    else:
        return class_to_int(code)

# The following classes support nicely loading and saving
# to disk in a format readable by postgres

class PGType(lazy_attribute):
    @classmethod
    def _load(cls, x):
        """
        Wraps :meth:`load` for the appropriate handling of NULLs.
        """
        if x != r'\N':
            return cls.load(xc)
    @classmethod
    def load(cls, x):
        """
        Takes a string from a file and returns the appropriate
        Sage object.

        Should be inverse to :meth:`save`

        This default function can be overridden in subclasses.
        """
        if x.isdigit():
            return ZZ(x)
        elif x.startswith('{'):
            return sage_eval(x.replace('{','[').replace('}',']'))
        else:
            return x

    @classmethod
    def _save(cls, x):
        """
        Wraps :meth:`save` for the appropriate handling of NULLs.
        """
        if x is None:
            return r'\N'
        else:
            return cls.save(x)

    @classmethod
    def save(cls, x):
        """
        Takes a Sage object stored in this attribute
        and returns a string appropriate to write to a file
        that postgres can load.

        Should be inverse to :meth:`load`

        This default function can be overridden in subclasses.
        """
        if isinstance(x, list):
            return str(x).replace('[','{').replace(']','}').replace("'",'"')
        else:
            return str(x)

class pg_text(PGType):
    pg_type = 'text'
class pg_smallint(PGType):
    pg_type = 'smallint'
class pg_integer(PGType):
    pg_type = 'integer'
class pg_numeric(PGType):
    # Currently only used to store large integers, so no decimal handling needed
    pg_type = 'numeric'
class pg_smallint_list(PGType):
    pg_type = 'smallint[]'
class pg_integer_list(PGType):
    pg_type = 'integer[]'
class pg_float8_list(PGType):
    pg_type = 'float8[]'
class pg_text_list(PGType):
    pg_type = 'text[]'
class pg_numeric_list(PGType):
    pg_type = 'numeric[]'
class pg_boolean(PGType):
    pg_type = 'boolean'
    @classmethod
    def load(cls, x):
        if x == 't':
            return True
        elif x == 'f':
            return False
        else:
            raise RuntimeError
    @classmethod
    def save(cls, x):
        if x:
            return 't'
        else:
            return 'f'
class pg_rational_list(PGType):
    pg_type = 'text[]'
    @classmethod
    def load(cls, x):
        def recursive_QQ(y):
            if isinstance(y, basestring):
                return QQ(y)
            else:
                return map(recursive_QQ, y)
        x = PGType.load(x)
        return recursive_QQ(x)
    @classmethod
    def save(cls, x):
        def recursive_str(y):
            if isinstance(y, list):
                return [recursive_str(z) for z in y]
            else:
                return str(y)
        x = recursive_str(x)
        return PGType.save(x)
class pg_jsonb(PGType):
    pg_type = 'jsonb'
    @classmethod
    def load(cls, x):
        return sage_eval(x)
    @classmethod
    def save(cls, x):
        return str(x).replace("'",'"')

class Stage(object):
    def __init__(self, controller, input, output):
        self.controller = controller
        self.input = input
        self.output = output

    @lazy_attribute
    def tasks(self):
        return [self.Task(g, q, self) for g, q in self.controler.gq]

class GenericTask(object):
    def __init__(self, g, q, stage):
        self.g, self.q, self.stage = g, q, stage
        self.logheader = stage.controller.logheader.format(g=g, q=q, name=stage.shortname)
    def ready(self):
        return all(ope(data[-1]) for data in self.input_data)
    def done(self):
        return all(ope(output.format(g=self.g, q=self.q)) for output in self.stage.output)
    @lazy_attribute
    def input_data(self):
        """
        Should be overridden in subclass
        """
        return []

class Worker(object):
    def __init__(self, logfile):
        self.logfile = logfile

class IsogenyClasses(object):
    """
    This class controls the creation of isogeny classes, grouped by g and q,
    as well as the process of loading and saving them to disk and to the LMFDB

    INPUT:

    - ``workers`` -- the number of processes to be allocated to this computation
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
        if config is None:
            if os.path.exists('config.ini'):
                config = os.path.abspath('config.ini')
            else:
                raise ValueError("Must have config.ini in directory or specify location")
        self.config_file = config
        cfgp = ConfigParser()
        cfgp.read(config)
        gs = sorted([ZZ(gx[1:]) for gx in config.options('extent')])
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

        # Create subdirectories if they do not exist
        basedir = os.path.abspath(os.path.expanduser(cfgp.get('dirs', 'base')))
        if not ope(basedir):
            os.makedirs(basedir)
        subdirs = [sub.strip() for sub in cfgp.get('dirs', 'subdirs').split(',')] + ['logs']
        for subdir in subdirs:
            if not ope(opj(basedir, subdir)):
                os.mkdir(opj(basedir, subdir))

        # Create Stages
        stages = []
        self.logfrequency = int(cfgp.get('logging', 'logfrequency'))
        self.logheader = cfgp.get('logging', 'logheader') + ' '
        for stage in cfgp.sections():
            if not stage.startswith('Stage'):
                continue
            info = [(key, cfgp.get(stage, key)) for key in cfgp.options(stage)]
            input = [val for (key, val) in info if key.startswith('in')]
            if any(key.startswith('out') for key, val in info):
                output, output_indexes = zip(*((val, key[3:]) for (key, val) in info if key.startswith('out')))
            else:
                output_indexes = output = []
            if any(key.startswith('data') for key, val in info):
                data, data_indexes = zip(*((val, key[4:]) for (key, val) in info if key.startswith('data')))
            else:
                data_indexes = data = []
            if output_indexes != data_indexes:
                raise ValueError("Output and data specifications need to be in the same order for %s.\nOne was %s while the other was %s" % (stage, ', '.join(output_indexes), ', '.join(data_indexes)))
            output = zip(output, data)
            stages.append(self.__class__.getattr(stage)(self, input=input, output=output))

        self.stages = stages
        self.tasks = sum((stage.tasks for stage in stages), [])
        logfile = cfgp.get('logging', 'logfile')
        self.workers = [Worker(opj(basedir, logfile.format(i=i))) for i in range(worker_count)]

    def run_serial(self):
        worker = self.workers[0]
        for task in self.tasks:
            task.run(worker.logfile)

    class StageGenerateSimple(Stage):
        name = 'Generate Simple'
        shortname = 'GenSimp'
        class Task(GenericTask):
            def run(self, logfile):
                stage = self.stage
                controller = stage.controller
                g, q = self.g, self.q
                with open(logfile, 'a') as logout:
                    def make_simples():
                        for i, Lpoly in enumerate(WeilPolynomials(2*g, q)):
                            IC = IsogenyClass(Lpoly=Lpoly)
                            try:
                                invs, mult = IC.simplepow_brauer_data
                            except ValueError:
                                continue
                            if mult != 1:
                                continue
                            yield IC
                    filename, attributes = self.stage.output
                    filename = filename.format(g=self.g, q=self.q)
                    controller.save(filename, make_simples(), attributes)

    class StageGenerateAll(Stage):
        name = 'Generate All'
        shortname = 'Gen All'
        class Task(GenericTask):
            @lazy_attribute
            def input_data(self):
                gs = range(1, self.g+1)
                q = self.q
                return [(g, self.input[0].format(g=g, q=q)) for g in gs]
            def run(self):
                stage = self.stage
                controller = stage.controller
                simples = {g: list(controller.load(filename)) for (g, filename) in self.input_data}
                def make_all():
                    for split in Partitions(g):
                        split_mD = multiplicity_dict(split)
                        it = cartesian_product_iterator([CombosWithReplacement(simples[g], c) for g, c in sorted(split_mD.items())])
                        for factors in it:
                            if len(factors) == 1:
                                factors = factors[0]
                            else:
                                factors = sum(factors, ())
                            factors = multiplicity_dict(factors).items()
                            yield IsogenyClass.by_decomposition(factors)
                filename, attributes = self.stage.output
                filename = filename.format(g=self.g, q=self.q)
                controller.save(filename, make_all(), attributes)
    class StageBasechange(Stage):
        name = 'Basechange'
        shortname = 'Bchange'
        @lazy_attribute
        def tasks(self):
            return [self.Task(g, q, self) for g, q in self.controler.gq if q.is_prime()]

        class Task(GenericTask):
            def __init__(self, g, p, stage):
                self.g, self.p, self.stage = g, p, stage
                controller = stage.controller
                self.qs = [q for (g, q) in controller.gq if q%p == 0]
                self.rs = [q.is_prime_power(get_data=True)[1] for q in self.qs]

            @lazy_attribute
            def input_data(self):
                return [(r, self.input[0].format(g=self.g, q=q)) for q, r in zip(self.qs, self.rs)]

            @lazy_attribute
            def output_data(self):
                stage = self.stage
                outputs = []
                for i, (filename, attributes) in enumerate(stage.output):
                    if i < 2:
                        # We use r=0 and r=-1 to encode the output files for endomorphism alg data
                        outputs.append((-i, filename.format(g=self.g, p=self.p), attributes))
                    else:
                        for r in self.rs:
                            outputs.append((r, filename.format(g=self.g, q=self.p^r), attributes))
                return outputs

            def run(self):
                stage = self.stage
                controller = stage.controller
                g, p = self.g, self.p
                in_db = {r: {IC.label: IC for IC in controller.load(filename)} for (r, filename) in self.input_data}
                geometric_degrees = default_dict(set)
                pair_lcms = defaultdict(set)
                for r, ICs in in_db.items():
                    for IC in ICs.values():
                        geometric_degrees[r].add(IC.geometric_extension_degree)
                    # Need to base change to all degrees in the database,
                    # even if these aren't pair lcms
                    for s in in_db.keys():
                        if s > r and s % r == 0:
                            pair_lcms[r].add(s // r)
                    for a, b in Combinations(geometric_degrees[r], 2):
                        pair_lcms[r].add(lcm(a,b))
                base_changes = defaultdict(lambda: defaultdict(list))
                simple_factors = []
                multiplicity_records = []
                factordata = [multiplicity_records, simple_factors]
                def compute_endomorphism_algebra(extension_class, base_class):
                    for simple_factor, mult in zip(extension_class.simple_distinct, extension_class.simple_multiplicities):
                        BCR = BaseChangeRecord(base_class, simple_factor, mult)
                        multiplicity_records.append(BCR)
                        simple_factors.append(simple_factor)
                for r, ICs in in_db.items():
                    lcms = pair_lcms[r]
                    q = p^r
                    for IC in ICs.values():
                        compute_endomorphism_algebra(IC)
                        geom_degree = IC.geometric_extension_degree
                        for s in lcms:
                            BC = IsogenyClass(Lpoly=base_change(IC.Lpoly, s, g=g, q=q))
                            base_changes[BC][r].append(IC)
                            if geom_degree % s == 0:
                                compute_endomorphism_algebra(BC)
                                if s == geom_degree:
                                    # record some geometric data
                                    IC.geometric_center_dim = BC.center_dim
                for r, ICs in in_db.items():
                    if r == 1:
                        continue
                    for IC in ICs.values():
                        if IC in base_changes:
                            models = sum(base_changes[IC].values(), [])
                            IC.primitive_models = []
                            for model in models:
                                if model not in base_changes:
                                    IC.primitive_models.append(model.label)
                        else:
                            IC.primitive_models = None
                        IC.twists = [] # filled in below
                for BC, ICs in base_changes.items():
                    for r in self.rs:
                        if len(ICs[r]) > 1: # rules out r=1
                            for IC in ICs[r]:
                                for JC in ICs[r]:
                                    if IC == JC:
                                        continue
                                    IC.twists.append(JC.label, BC.label, BC.r // r)
                for r, filename, attributes in self.output_data:
                    if r <= 0:
                        controller.save(filename, factordata[r], attributes)
                    else:
                        controller.save(filename, in_db[r], attributes)
    class StageCombine(Stage):
        name = 'Combine'
        shortname = 'Combine'
        class Task(GenericTask):
            @lazy_attribute
            def input_data(self):
                return [filename.format(g=self.g, q=self.q) for filename in self.stage.input]
            def run(self):
                stage = self.stage
                controller = stage.controller
                sources = [controller.load(filename) for filename in self.input_data]
                data = defaultdict(list)
                for source in sources:
                    for IC in source:
                        data[IC.label].append(IC)
                def make_all():
                    for key in sorted(data.keys()):
                        ICs = data[key]
                        yield IsogenyClass.combine(*ICs)
                outfile, attributes = self.output[0]
                outfile = outfile.format(g=self.g, q=self.q)
                controller.save(outfile, make_all(), attributes)
    @staticmethod
    def load(filename, start=None, stop=None, cls=None):
        """
        Iterates over all of the isogeny classes stored in a file.
        The data contained in the file is specified by header lines: the first giving the
        column names (which are identical to PGType attributes of IsogenyClass)
        the second giving postgres types (unused here), and the third blank.

        We use ':' as a separator.

        INPUT:

        - ``filename`` -- the filename to open
        - ``start`` -- if not ``None``, the line to start from, indexed so that 0 is
            the first line of data.  Note that this will be the fourth line of the
            file because of the header.
        - ``stop`` -- if not None``, the first line not to read, indexed as for ``start``.
        """
        if cls is None:
            cls = IsogenyClass
        with open(filename) as F:
            for i, line in enumerate(F):
                if i == 0:
                    header = line.strip().split(':')
                elif i >= 3 and (start is None or i-3 >= start) and (stop is None or i-3 < start):
                    yield IsogenyClass.load(line.strip(), header)

    @staticmethod
    def save(filename, isogeny_classes, attributes, cls=None, force=False):
        """
        INPUT:

        - ``filename`` -- a filename to write
        - ``isogeny_classes`` -- an iterable of instances to write to the file
        - ``attributes`` -- a list of attributes to save to the file
        - ``cls`` -- the class of the entries of ``isogeny_classes``
        - ``force`` -- if True, will allow overwriting an existing file
        """
        if not force and ope(filename):
            raise ValueError("File %s already exists")
        if cls is None:
            cls = IsogenyClass
        types = [getattr(cls, attr) for attr in attributes]
        header = [':'.join(attributes),
                  ':'.join(attr.pg_type for attr in types),
                  '\n']
        with open(filename, 'w') as F:
            F.write('\n'.join(header))
            for isog in isogeny_classes:
                F.write(isog.save(attributes) + '\n')

    # make sure to sort results appropriately

class PGSaver(object):
    @classmethod
    def load(cls, s, header):
        """
        INPUT:

        - ``s`` -- a string, giving the data defined by ``header`` as colon separated values
        - ``header`` -- a list of attribute names to fill in
        """
        data = s.split(':')
        isoclass = cls()
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
        return ':'.join(getattr(cls, attr)._save(getattr(self, attr)))

class IsogenyClass(PGSaver):
    """
    An isogeny class of abelian varieties over a finite field, as constructed from a Weil polynomial
    """
    @classmethod
    def by_label(cls, label, check=True):
        if label.count('.') != 2:
            raise ValueError("Invalid label")
        g, q, wp = label.split('.')
        try:
            g = int(g)
            q = int(q)
        except ValueError:
            raise ValueError("Invalid label")
        if wp.count('_') != g-1:
            raise ValueError("Invalid label")
        coeffs = map(signed_class_to_int, wp.split('_'))
        coeffs = [1r] + coeffs + [q^i*c for i,c in enumerate(reversed(coeffs[:-1]))] + [q^g]
        Lpoly = ZZ['x'](coeffs)
        return cls(Lpoly, check)

    @classmethod
    def by_decomposition(cls, factors):
        """
        INPUT:

        - ``factors`` -- a list of pairs (IC, e) where IC is a simple isogeny class and e is an exponent
        """
        Lpoly = prod(IC.Lpoly^e for IC,e in factors)
        result = cls(Lpoly=Lpoly)
        result.decomposition = factors
        result.has_decomposition = True
        return result

    def __init__(self, Lpoly=None, poly=None, label=None):
        # All None is allowed since the load method writes to the fields outside the __init__ method
        if Lpoly is not None:
            self.Lpoly = Lpoly
        if poly is not None:
            self.poly = poly
        if label is not None:
            self.label = label
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
        coeffs = [1] + coeffs + [q^i * c for (i, c) in enumerate(reversed(coeffs[:-1]))] + [q^g]
        return coeffs

    @lazy_attribute
    def Lpoly(self):
        return ZZ['x'](self.poly)

    @pg_text
    def label(self):
        return '%s.%s.%s' % (
            self.g,
            self.q,
            '_'.join(signed_cremona_letter_code(c) for c in self.Lpoly[1:g+1]))

    @lazy_attribute
    def Ppoly(self):
        return self.Lpoly.reverse()

    @pg_smallint
    def g(self):
        d = self.Lpoly.degree()
        if d % 2 == 1:
            raise ValueError("Must have even degree")
        return d // 2

    @pg_integer
    def q(self):
        return self.Ppoly[0].nth_root(self.g)

    @lazy_attribute
    def p(self):
        p, r = q.is_prime_power(get_data=True)
        self.r = r
        return p

    @lazy_attribute
    def r(self):
        p, r = q.is_prime_power(get_data=True)
        self.p = p
        return r

    @lazy_attribute
    def Ppoly_factors(self):
        return self.Ppoly.factor()

    @pg_rational_list
    def slopes(self):
        p, r, g = self.p, self.r, self.g
        np = self.Lpoly.change_ring(Qp(p)).newton_polygon()
        return [(np(i) - np(i+1)) / r for i in range(2*g)]

    @pg_smallint
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
        return self._real_poly(self.Lpoly, self.g, self.q)

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
        return f.log().derivative().coefficients()[:prec]

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
    def point_count(self):
        return self.curve_counts[0]

    @pg_text
    def poly_str(self):
        return ' '.join(map(str, self.poly))

    @pg_text
    def abvar_counts_str(self):
        return ' '.join(map(str, self.abvar_counts))

    @pg_text
    def curve_counts_str(self):
        return ' '.join(map(str, self.curve_counts))

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
        if not self.is_simple:
            raise ValueError("Non-simple")
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
        a = self.K.gen()
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
            return Factorization([(IsogenyClass(label=label), e) for label, e in zip(self.simple_distinct, self.simple_multiplicities)])
        else:
            # Lazy attributes store their result as an attribute, so the only way we can reach this branch
            # is if the other recursed back.  So we compute the decomposition from the factorization.
            factorization = []
            for factor, power in self.Ppoly_factors:
                IC = IsogenyClass(factor^power)
                invs, e = IC.simplepow_brauer_data
                if e != 1:
                    IC = IsogenyClass(factor^(power//e))
                    IC.is_simple = True
                factorization.append((IC, e))
            factorization.sort(key=lambda pair: (pair[0].g, pair[0].poly))
            return Factorization(factorization, sort=False)

    @pg_boolean
    def is_simple(self):
        return len(self.decomposition) == 1 and self.decomposition[0][1] == 1

    @pg_text_list
    def simple_distinct(self):
        return [IC.label for IC, e in self.decomposition]

    @pg_smallint_list
    def simple_multiplicities(self):
        return [e for IC, e in self.decomposition]

    @pg_rational_list
    def brauer_invariants(self):
        if not self.is_simple:
            raise ValueError("Non-simple")
        return self.simplepow_brauer_data[0]

    @pg_rational_list
    def places(self):
        """
        A list of rational numbers, giving coordinates for the second entry
        of the two-generator representation of the ideals above p.

        Only used for simple isogeny classes
        """
        if not self.is_simple:
            raise ValueError("Non-simple")
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
        for poly, e in self.Ppoly_factors:
            coeffs = R(pari(poly).polredbest().polredabs()).coefficients(sparse=False)
            rec = db.nf_fields.lucky({'coeffs':coeffs}, projection=['label','degree','galt'], sort=[])
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
        return self._nf_data[0][0]

    @pg_text
    def center(self):
        # Used as a synonym when storing endomorphism data for simple factors of base changes
        if not self.is_simple:
            raise ValueError("Non-simple")
        return self.number_field

    @pg_smallint
    def center_dim(self):
        if self.has_decomposition:
            return sum(simple_factor.center_dim for simple_factor, e in self.decomposition)
        else:
            return sum(poly.degree() for poly, e in self.Ppoly_factors)

    @pg_smallint
    def geometric_center_dim(self):
        g, q, s = self.g, self.q, self.geometric_extension_degree
        return IsogenyClass(Lpoly=base_change(self.Lpoly, s, g=g, q=q)).center_dim

    @pg_text
    def galois_group(self):
        # Used for caching after the first stage
        if not self.is_simple:
            raise ValueError("Non-simple")
        return self._nf_data[1][0]

    @pg_text
    def number_fields(self):
        return self._nf_data[0]

    @pg_text
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

    @lazy_attribute
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

        fieldext = ZZ(1)
        for factor, power in square.factor():
            m = Cyclotomic().order(factor)
            if m > 0:
                fieldext = fieldext.lcm(m)

        return fieldext

    @pg_text_list
    def primitive_models(self):
        # Will be set by the base change stage code
        raise RuntimeError

    @pg_boolean
    def is_primitive(self):
        return self.primitive_models is None

    @pg_jsonb
    def twists(self):
        # Will be set by the base change stage code
        raise RuntimeError


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
        d, a, b = xgcd(f.change_ring(QQ), g.change_ring(QQ))
        if d.degree() > 0:
            return 0
        n = lcm(c.denominator() for c in a)
        h = h1*h2
        H, rem = h.quo_rem(x^2-4*q)
        if rem == 0:
            splitelt = n*a1*h1 # 0 mod h1, n mod h2
            otherelt = splitelt + x*H
            g = gcd(ZZ(c) for c in otherelt)
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
        elif g % 2 == 1:
            # Every odd-dimensional simple isogeny class has a principal polarization
            return 1r
        elif self.is_simple:
            plus_poly = self.real_Lpoly
            # Look at the CM field K and its real subfield K+.
            # If K/K+ is ramified at a finite prime, or if there is a prime of K+
            # that divides F - V and that is inert in K/K+,
            # then there's a principally polarized variety in the isogeny class.
            poly = plus_poly.parent()(coeffs)
            K.<F> = NumberField(poly)
            Kplus.<Fplus> = NumberField(plus_poly)
            D = K.discriminant()
            Dplus = Kplus.discriminant()
            if D.abs() != Dplus^2:
                return 1r
            V = q / F
            # F -> V is complex conjugation
            conj = K.hom([V])
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
        elif all(IC.has_principal_polarization for IC in self.decomposition):
            # If every factor can be principally polarized, then so can the product
            # The converse isn't true
            return 1r
        return 0r

    @pg_smallint
    def has_jacobian(self):
        g, q, p, r = self.g, self.q, self.p, self.r
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
        this isogeny class cannot contain a Jacobian.
        """
        return any(cnt < 0 for cnt in self._newpoints)

    # The following _nojac methods are from Howe's http://ewhowe.com/Magma/IsogenyClasses.magma

    def _nojac_stohr_voloch(self):
        # For certain values of q and g, the paper [Stohr and Voloch 1986]
        # gives bounds on N_q(g) that may be better than the Oesterle bound.
        g, q, p = self.g, self.q, self.p
        if g < 3 or p < 2*g - 3:
            # No info in this case
            return
        N = self.point_count

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
                res = self.modified_reduced_resultant(f, g)
                if res == 1:
                    return True

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

for d in range(1,6):
    def dim_factors(self):
        return sum([e for (simp_label, e) in zip(self.simple_distinct, self.simple_multiplicities)
                    if simp_label.split('.')[0] == str(d)])
    def dim_distinct(self):
        return len([simp_label for simp_label in self.simple_distinct if simp_label.split('.')[0] == str(d)])
    setattr(IsogenyClass, 'dim%d_factors'%d, pg_smallint(dim_factors))
    setattr(IsogenyClass, 'dim%d_distinct'%d, pg_smallint(dim_distinct))

class BaseChangeRecord(object):
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

def angle_rank(u,p):
    """
    There are two methods for computing this.
    The first method is to use S-units in sage where S = primes in Q(pi) dividing p.
    The second method is to use lindep in pari.
    """

    """
    K.<a> = u.splitting_field()
    l = [p] + [i[0] for i in u.roots(K)]
    S = K.primes_above(p)
    UGS = UnitGroup(K, S = tuple(S), proof=False)
    ## Even with proof=False, it is guaranteed to obtain independent S-units; just maybe not the fully saturated group.
    d = K.number_of_roots_of_unity()
    gs = [K(i) for i in UGS.gens()]
    l2 = [UGS(i^d).exponents() for i in l] #For x = a^1b^2c^3 exponents are (1,2,3)
    for i in range(len(l)):
        assert(l[i]^d == prod(gs[j]^l2[i][j] for j in range(len(l2[i]))))
    M = Matrix(l2)
    return M.rank()-1
    """

oldmatcher = re.compile(r"weil-(\d+)-(\d+)\.txt")
simplematcher = re.compile(r"weil-simple-g(\d+)-q(\d+)\.txt")
allmatcher = re.compile(r"weil-all-g(\d+)-q(\d+)\.txt")
LoadedPolyData = namedtuple("LoadedPolyData","label poly angle_numbers p_rank slopes invs places")
def load_previous_polys(q = None, g = None, rootdir=None, all = False):
#needs the fields to have the proper quotations/other Json formatting.
    if rootdir is None:
        rootdir = os.path.abspath(os.curdir)
    if all:
        matcher = allmatcher
    else:
        matcher = simplematcher
    D = defaultdict(list)
    R = PolynomialRing(QQ,'x')
    def update_dict(D, filename):
        with open(filename) as F:
            for line in F.readlines():
                data = json.loads(line)
                label, g, q, polynomial, angle_numbers = data[:5]
                p_rank, slopes = data[6:8]
                invs, places = data[13:15]
                D[g,q].append(LoadedPolyData(label, R(polynomial), angle_numbers, p_rank, slopes, invs, places))
    if q is not None and g is not None:
        filename = "weil-simple-g%s-q%s.txt"%(g, q)
        update_dict(D, filename)
    else:
        for filename in os.listdir(rootdir):
            match = matcher.match(filename)
            if match:
                gf, qf = map(int,match.groups())
                if q is not None and qf != q:
                    continue
                if g is not None and gf != g:
                    continue
                update_dict(D, os.path.join(rootdir, filename))
    return D

def alternating(pol, m):
    """
    This appears to take forever but there is a precomputed version elsewhere.
    """
    d = pol.degree()
    pl = pol.list()
    e = SymmetricFunctions(QQ).e()
    dm = binomial(d, m)
    l = [(-1)^i*e[i](e[m]).restrict_parts(d) for i in range(dm+1)]
    P = pol.parent()
    R = P.base_ring()
    ans = []
    for i in range(dm,-1,-1):
        s = R.zero()
        u = tuple(l[i])
        for j, c in u:
            s += R(c) * prod(pl[d-k] for k in j)
        ans.append(s)
    return P(ans)

def find_invs_and_slopes(p,r,P):
    ### KEEPS OLD INVARIANTS BEHAVIOR (doesn't reduce mod Z) ###
    poly = P.change_ring(QQ)
    K.<a> = NumberField(poly)
    l = K.primes_above(p)
    invs = []
    slopes = []
    for v in l:
        vslope = a.valuation(v)/K(p^r).valuation(v)
        slopes.append(vslope)
        vdeg = v.residue_class_degree()*v.ramification_index()
        invs.append(vslope*vdeg)
    return invs,slopes

