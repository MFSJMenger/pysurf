from numpy import (isscalar, zeros, zeros_like, dot, 
                   array, ndarray, complex128, outer, diag, linalg, exp, 
                   abs, imag, real, maximum, sqrt, cumsum, less, copy, sum, pi)
from random import uniform
from time import (time, ctime)
from collections import namedtuple
from abc import abstractmethod
from pysurf.spp import SurfacePointProvider
from pysurf.database import PySurfDB
from colt import Colt

class VelocityVerletPropagator:

    def __init__(self, state):
        self.state = state
        self.t = self.state.t
        self.dt = self.state.dt
        self.t_max = self.dt*self.state.mdsteps 
        if self.state.method == "Surface_Hopping":
            self.electronic = SurfaceHopping(self.state)
            if self.state.ncoeff[self.state.instate] == 0:
                raise SystemExit("Wrong population for initial state")
        elif self.state.method == "Born_Oppenheimer":  
            self.electronic = BornOppenheimer(self.state)
        self.results = PrintResults()

    def run(self):
        
        if self.t > self.t_max:
            raise SystemExit("Noting to be done")
        
        state = self.state       
        results = self.results 
        grad_old = self.electronic.setup(state)
        acce_old = self.accelerations(state, grad_old)

        results.print_head(state)
        while(self.t <= self.t_max):
            """updating coordinates"""
            crd_new = self.positions(state, acce_old, self.dt)
            """updating accelerations"""    
            grad_new = self.electronic.new_surface(state, results, crd_new, self.t, self.dt)
            acce_new = self.accelerations(state, grad_new)
            """updating velocities"""
            vel_new = self.velocities(state, acce_old, acce_new, self.dt)
            """updating variables"""
            acce_old = self.update_state(state, acce_new, crd_new, vel_new) 
            self.t += self.dt 
        results.print_bottom(state)

    def accelerations(self, state, grad):
        if isscalar(state.mass) and isscalar(grad[state.instate]):
            return -grad[state.instate]/state.mass
        else:
            gradient = grad[state.instate]
            acce = zeros(gradient.shape)
            for i,m in enumerate(state.mass):
                acce[i] = -gradient[i]/m
            return acce

    def positions(self, state, a_0, dt):
        return state.crd + state.vel*dt + 0.5*a_0*dt**2

    def velocities(self, state, a_0, a_1, dt):
        return state.vel + 0.5*(a_0 + a_1)*dt 

    def update_state(self, state, acce_new, crd_new, vel_new):
        state.crd = crd_new
        state.vel = vel_new
        state.e_two_prev_steps = state.e_prev_step
        state.e_prev_step = state.e_curr
        state.e_curr = state.ene
        acce_old = acce_new
        return acce_old

class BornOppenheimer:

    needed_properties = ["energy","gradient"]

    def __init__(self, state):
        self.nstates = state.nstates
        self.natoms = state.natoms
        self.spp = SurfacePointProvider.from_questions(["energy","gradient"], self.nstates, self.natoms, config ="spp.inp", atomids = state.atomids, check_only=True)

    @abstractmethod
    def get_gradient(self, crd, curr_state):
        result = self.spp.request(crd, ['gradient'], states=[curr_state])
        return result['gradient']

    def get_energy(self, crd):
        result = self.spp.request(crd, ['energy'])
        return result['energy']

    def cal_ekin(self, mass, vel):
        ekin = 0
        if isscalar(mass) and isscalar(vel):
            ekin = 0.5*mass*vel**2
        else:
            for i, m in enumerate(mass):
                ekin += 0.5*m*dot(vel[i],vel[i])
        return ekin
    
    def setup(self, state):
        state.ene = self.get_energy(state.crd)
        grad = self.get_gradient(state.crd, state.instate)
        state.epot = state.ene[state.instate]
        state.ekin = self.cal_ekin(state.mass, state.vel)
        return grad

    def new_surface(self, state, results, crd_new, t, dt):
        grad_new = self.get_gradient(crd_new, state.instate)
        results.print_bh_var(t, dt, state) #printing variables 
        results.save_db(t,state) #save variables in database
        state.ene = self.get_energy(crd_new)
        state.epot = state.ene[state.instate]
        state.ekin = self.cal_ekin(state.mass, state.vel)
        return grad_new


class Propagator:

    def __init__(self, state):
        self.prob_name = state.prob
        self.nstates = state.nstates

    def elec_density(self, state):
        c_mch = state.ncoeff
        if isinstance(c_mch, ndarray) != True:
            c_mch = array(c_mch,dtype=complex128) 
        return outer(c_mch, c_mch.conj())

    def grad(self, ene_cou_grad):
        g_mch = ene_cou_grad.grad
        if self.prob_name == "diagonal":
            u = ene_cou_grad.u
            g_diag = {} #it must request all the gradients for all the electronic states in diagonal approach 
            for i in range(self.nstates):
                g_diag.update({i:dot(u.T.conj()[i,:],u[:,i]).real*g_mch[i]})
            return g_diag
        elif self.prob_name in  ("tully", "lz"):
            return g_mch

    def mch_propagator(self, h_mch, vk, dt):
        h_total = diag(h_mch) - 1j*(vk) 
        ene, u = linalg.eigh(h_total)
        p_mch = linalg.multi_dot([u, diag(exp( -1j * ene * dt)), u.T.conj()])
        return p_mch

    def elec_density_new(self, state,  rho_old, p_mch):
        """ Computing propagation and new density:
            D and U are diagonal and unitary matrices of hamiltonian 
            rho_ij(t) = c_i(t)*c_j(t)^{*}
            as c(t+dt) = U*exp(-j*D*dt)*U.T*c(t),
            then rho_ij(t+dt) = c_i(t+dt)*c_j(t+dt)^{*}
                              = U*exp(-j*D*dt)*U.T*c_i(t)*c_j(t)^{*}*U*exp(j*D*dt)*U.T
                              = U*exp(-j*D*dt)*U.T*rho_ij(t)*U*exp(j*D*dt)*U.T
        """
        c_dt = linalg.multi_dot([p_mch, state.ncoeff]) 
        rho_new = linalg.multi_dot([p_mch, rho_old, p_mch.T.conj()])   
        c_rho_new = namedtuple("c_rho_new", "c_dt rho_new") 
        return c_rho_new(c_dt, rho_new)

    def hopping_probability(self, c_j_dt, c_i_dt, c_i_t, p_diag_ji, p_diag_ii): 
        prob_factor_1 = 1 - abs(dot(c_i_dt, c_i_dt.conj()))/abs(dot(c_i_t, c_i_t.conj())) 
        prob_factor_2_N = ( dot(c_j_dt, dot(p_diag_ji.conj(), c_i_t.conj()) ) ).real
        prob_factor_2_D = ( abs(dot(c_i_t, c_i_t.conj())) -\
                       ( dot(c_i_dt, dot(p_diag_ii.conj(), c_i_t.conj()) ) ).real)
        if prob_factor_2_D == 0:
            prob_ji = 0.0
        elif prob_factor_1*((prob_factor_2_N/prob_factor_2_D)) < 0:
            prob_ji = 0.0
        else:
            prob_ji = prob_factor_1*((prob_factor_2_N/prob_factor_2_D))
        return prob_ji

    def probabilities_tully(self, state, vk_old, ene_new, dt):
        rho_old = self.elec_density(state)
        instate = state.instate
        ene_old = state.ene
        vk_new = state.vk
        ene = 0.5*(ene_old + ene_new)
        vk = 0.5*(vk_old + vk_new)
        h_total = diag(ene) - 1j*(vk)
        p_mch = self.mch_propagator(ene, vk, dt)
        probs = (2.0 * imag(rho_old[instate,:] * h_total[:,instate]) * dt)/(real(rho_old[instate,instate])) 
        probs[instate] = 0.0
        probs = maximum(probs, 0.0)
        tully = namedtuple("tully", "probs rho_old p_mch")
        return tully(probs, rho_old, p_mch) 
 
    def probabilities_diagonal(self, state, diag_prop):
        c_diag_dt = diag_prop.c_diag_new
        c_diag = diag_prop.c_diag
        p_diag_dt = diag_prop.p_diag_new
        probs = zeros(self.nstates)
        instate = state.instate
        for i in range(self.nstates):
            probs[i] = self.hopping_probability(c_diag_dt[i], c_diag_dt[instate],\
                                                c_diag[instate],p_diag_dt[i,instate], p_diag_dt[instate,instate])
        probs[instate] = 0.0
        return probs 

    def compute_diff(self,e,i,j):
        return e[i]-e[j]

    def abs_gt(self, a, b):
        return abs(a) > abs(b)

    def lz_select(self, state, dt):
        """Compute the hopping in LZ approach"""
        if state.e_two_prev_steps is None:
            iselected = None
            prob = -1.0
        else:
            e_curr = state.e_curr
            e_two_step_prev = state.e_two_prev_steps 
            e_prev_step = state.e_prev_step
            iactive = state.instate
            iselected = None
            prob = -1.0
            for istate in range(self.nstates):
                #
                if istate == iactive:
                    continue
                # compute energy differences
                d_two_prev_steps = self.compute_diff(e_two_step_prev, iactive, istate)
                d_prev_step = self.compute_diff(e_prev_step, iactive, istate)
                d_curr = self.compute_diff(e_curr, iactive, istate)
                # compute probability
                if (self.abs_gt(d_curr, d_prev_step) and self.abs_gt(d_two_prev_steps, d_prev_step)):
                    curr_hop_prob = self.probabilities_lz(d_curr, d_prev_step, d_two_prev_steps, dt)
                    if curr_hop_prob > prob:
                        prob = curr_hop_prob
                        iselected = istate
        ise_prop = namedtuple("ise_prop","iselected prob")
        return ise_prop(iselected, prob)

    def probabilities_lz(self, d_curr, d_prev_step, d_two_prev_steps, dt):
        """Compute the hopping probability between two electronic states"""
        finite_difference_grad = ((d_curr + d_two_prev_steps - 2 * d_prev_step) / dt**2)
        return exp((-pi/2.0) * sqrt(abs(d_prev_step)**3 / abs(finite_difference_grad)))

    def diag_propagator(self, ene_cou_grad, dt, state):
        c_mch = state.ncoeff
        u_new = ene_cou_grad.u
        u = state.u
        ene = state.ene
        vk = state.vk
        if isinstance(c_mch, ndarray) != True:
            c_mch = array(c_mch)
        c_diag = dot(u.T.conj(),c_mch)
        p_mch_new = self.mch_propagator(ene, vk, dt)
        p_diag_new = dot(u_new.T.conj() ,dot(p_mch_new,u))
        c_diag_new = dot(p_diag_new,c_diag) 
        diag_prop = namedtuple("diag_prop","c_diag, c_diag_new, p_diag_new")
        return diag_prop(c_diag, c_diag_new, p_diag_new)
 
    def new_prob_grad(self, state, vk_old, ene_cou_grad, dt):
        if self.prob_name == "diagonal":
            grad_new = self.grad(ene_cou_grad)
            diag_prop = self.diag_propagator(ene_cou_grad, dt, state)
            probs = self.probabilities_diagonal(state, diag_prop)
            result = namedtuple("result","probs grad_new diag_prop") 
            return result(probs, grad_new, diag_prop)
        elif self.prob_name == "tully":
            #tully = self.probabilities_tully(state, dt)
            ene_new = ene_cou_grad.ene
            tully = self.probabilities_tully(state, vk_old, ene_new, dt)
            probs = tully.probs
            grad_new = ene_cou_grad.grad 
            result = namedtuple("result","probs grad_new tully") 
            return result(probs, grad_new, tully)
        elif self.prob_name == "lz":
            ene_new = ene_cou_grad.ene
            lz = self.lz_select(state, dt)
            grad_new = ene_cou_grad.grad 
            result = namedtuple("result","grad_new lz") 
            return result(grad_new, lz)
        else:
            raise SystemExit("A right probability method is not defined")

    def new_hopp_grad(self, state, u_grad):
        if self.prob_name == "diagonal":
            grad_new = self.grad(u_grad)
            return grad_new 
        elif self.prob_name == "tully":
            grad_new = u_grad.grad 
            return grad_new 
        else:
            raise SystemExit("A right probability method is not defined")

    def instantaneous_decoherence_correction(self, state, ncoeff):
        if isinstance(ncoeff, ndarray) != True:
            ncoeff = array(ncoeff)
        ncoeff = zeros_like(ncoeff)
        ncoeff[state.instate] = 1
        state.ncoeff = ncoeff

    def energy_based_decoherence_correction(self, state, ncoeff):
        const = 1.0 + (0.1/state.ekin)
        curr = state.instate
        add = 0.0
        ncoeff_prima = [0.0]*state.nstates
        for k in range(state.nstates):
            if k != curr:
                t_kl = const/abs(state.ene[curr]-state.ene[k])
                ncoeff_prima[k] = ncoeff[k]*exp(-0.5*state.dt/t_kl)
                add += abs(dot(ncoeff_prima[k], ncoeff_prima[k].conj()))
        ncoeff_prima[curr] = ncoeff[curr]*sqrt((1.0 - add)/(abs(dot(ncoeff[curr], ncoeff[curr].conj()))))
        if isinstance(ncoeff_prima, ndarray) != True:
            ncoeff_prima = array(ncoeff_prima)
        state.ncoeff = ncoeff_prima 

    def no_decoherence_correction(self, state, ncoeff):
        if isinstance(ncoeff, ndarray) != True:
            ncoeff = array(ncoeff)
        state.ncoeff = ncoeff 

    def check_coherence(self, state, ncoeff, hop, att, succ):
        if state.decoherence == "No_DC":
            self.no_decoherence_correction(state, ncoeff)
        elif state.decoherence == "EDC":
            self.energy_based_decoherence_correction(state, ncoeff)
        elif state.decoherence == "IDC_A":
            if att =="yes":
                self.instantaneous_decoherence_correction(state, ncoeff)
            else:
                self.no_decoherence_correction(state, ncoeff)
        elif state.decoherence == "IDC_S":
            if succ =="yes":
                self.instantaneous_decoherence_correction(state, ncoeff)
            else:
                self.no_decoherence_correction(state, ncoeff)
        else:
            raise SystemExit("A right decoherence method is not defined")  

    def new_ncoeff(self, state, grad_probs, hop, att, succ):
        if self.prob_name == "diagonal":   
            ncoeff = dot(state.u, grad_probs.diag_prop.c_diag_new)
            self.check_coherence(state, ncoeff, hop, att, succ)
            state.rho = self.elec_density(state)
        elif self.prob_name == "tully":
            c_rho_new = self.elec_density_new(state, grad_probs.tully.rho_old, grad_probs.tully.p_mch)
            state.rho = c_rho_new.rho_new
            self.check_coherence(state, c_rho_new.c_dt, hop, att, succ)
        else:
            raise SystemExit("A right probability method is not defined")  

class RescaleVelocity:

    def __init__(self, state, ene_cou_grad):
        self.coupling = state.coupling
        self.rescale_vel = state.rescale_vel
        self.mass = state.mass
        self.state_old = state.instate
        self.ene_new = ene_cou_grad.ene
        if self.rescale_vel == "nacs" ==  self.coupling:
            self.nac_new = ene_cou_grad.nac
            self.nac_old = state.nac

    def direction(self, vel, state_new):
        if self.rescale_vel == "nacs" ==  self.coupling:
            return (0.5)*(self.nac_old[state_new,self.state_old] + self.nac_new[state_new,self.state_old])
        elif self.rescale_vel == "momentum":
            p = zeros(vel.shape)
            for i, m in enumerate(self.mass):
                p[i] = vel[i]*m
            return p
        
         
    def diff_ji(self, state_new):
        return self.ene_new[self.state_old] - self.ene_new[state_new]
    
    def beta_ji(self, vel, direct):
        if isscalar(vel):
            return vel*direct
        return dot(vel.flatten(),direct.flatten())
    
    def alpha_ji(self, direct):
        if isscalar(self.mass):
            return 0.5*(direct**2)/self.mass
        else:
            alpha = 0.0
            for i, m in enumerate(self.mass):
                alpha += dot(direct[i], direct[i])/m
            return 0.5*alpha

    def new_velocity(self, state, gama_ji, direct):
        if isscalar(state.vel) and isscalar(self.mass):
            state.vel = state.vel - gama_ji*(direct/self.mass)
        else:
            for i, m in enumerate(self.mass):
                state.vel[i] = state.vel[i] - gama_ji*(direct[i]/m)

    def rescale_velocity(self, state, state_new):
        direct = self.direction(state.vel, state_new)
        diff = self.diff_ji(state_new)
        beta = self.beta_ji(state.vel, direct)
        alpha = self.alpha_ji(direct) 
        if (beta**2 + 4*alpha*diff) < 0.0:
            """
            If this condition is satisfied, there is not hopping and 
            then the nuclear velocity simply are reversed.
            """
            gama_ji = beta/alpha
            self.new_velocity(state, gama_ji, direct)
            hop = "not"
            return hop 
        else:
            """
            If this condition is satisfied, a hopping from 
            current state to the first true state takes place 
            and the current nuclear velocity is ajusted in order 
            to preserve the total energy.
            """
            if beta < 0.0:
                gama_ji = (beta + sqrt(beta**2 + 4*alpha*diff))/(2*alpha)
                self.new_velocity(state, gama_ji, direct)
            else:
                gama_ji = (beta - sqrt(beta**2 + 4*alpha*diff))/(2*alpha)
                self.new_velocity(state, gama_ji, direct)
            hop = "yes"
            return hop


class SurfaceHopping(BornOppenheimer):

    def __init__(self, state):
        self.nstates = state.nstates
        self.mass = state.mass
        self.natoms = state.natoms
        self.coupling = state.coupling
        self.vel_old = zeros_like(state.vel) 
        self.prob = state.prob       
        if self.coupling == "nacs" and self.prob == "tully":
            needed_properties = ["energy", "gradient", "nacs"]
            self.spp = SurfacePointProvider.from_questions(["energy","gradient","nacs"], self.nstates, self.natoms, config ="spp.inp", atomids = state.atomids)
        elif self.coupling == "wf_overlap" and self.prob == "tully":
            needed_properties = ["energy", "gradient", "wf_overlap"]
            self.spp = SurfacePointProvider.from_questions(["energy","gradient","wf_overlap"], self.nstates, self.natoms, config ="spp.inp", atomids = state.atomids)
        elif self.coupling == "non_coup" and self.prob == "lz":
            needed_properties = ["energy", "gradient"]
            self.spp = SurfacePointProvider.from_questions(["energy","gradient"], self.nstates, self.natoms, config ="spp.inp", atomids = state.atomids)

    def get_gradient(self, crd, curr_state):
        result = self.spp.request(crd, ['gradient'], states=[curr_state])
        return result['gradient']

    def get_energy(self, crd):
        result = self.spp.request(crd, ['energy'])
        return result['energy']

    def get_coupling(self, crd):
        if self.coupling == "nacs":
            result = self.spp.request(crd, ['nacs'])
            return result['nacs']
        elif self.coupling == "wf_overlap":
            result = self.spp.request(crd, ['wf_overlap'])
            return result['wf_overlap']
            
    def get_ene_cou_grad(self, crd, curr_state):
        grad = self.get_gradient(crd, curr_state) 
        h_mch = self.get_energy(crd)
        ene, u = linalg.eigh(diag(h_mch))
        if self.coupling == "nacs":
            nac = self.get_coupling(crd)
            ene_cou_grad = namedtuple("ene_cou_grad", "ene u nac grad")
            return ene_cou_grad(ene, u, nac, grad)
        elif self.coupling == "wf_overlap":
            wf_ov = self.get_coupling(crd)
            ene_cou_grad = namedtuple("ene_cou_grad", "ene u wf_ov grad")
            return ene_cou_grad(ene, u, wf_ov, grad)
        elif self.coupling == "non_coup":
            ene_cou_grad = namedtuple("ene_cou_grad", "ene grad")
            return ene_cou_grad(ene, grad)

    def get_hopp_u_grad(self, crd, state):
        grad = self.get_gradient(crd, state.instate)
        u = state.u
        u_grad = namedtuple("u_grad", "u grad")
        return u_grad(u, grad) 

    def vk_coupl_matrix(self, nac, vel):
        vk = zeros((self.nstates,self.nstates))
        if isscalar(vel):
            for i in range(self.nstates):
                for j in range(self.nstates):
                    if i < j:
                        vk[i,j] = vel*nac[i,j]
                        vk[j,i] = -vk[i,j]
        else:
            for i in range(self.nstates):
                for j in range(self.nstates):
                    if i < j:
                        vk[i,j] = dot(vel.flatten(),nac[i,j].flatten())
                        vk[j,i] = -vk[i,j]
        return vk

    def cal_ekin(self, mass, vel):
        ekin = 0
        if isscalar(mass) and isscalar(vel):
            ekin = 0.5*mass*vel**2
        else:
            for i, m in enumerate(mass):
                ekin += 0.5*m*dot(vel[i],vel[i])
        return ekin

    def setup(self, state):
        ene_cou_grad = self.get_ene_cou_grad(state.crd, state.instate)
        propagator = Propagator(state)
        grad_old = propagator.grad(ene_cou_grad)
        state.ene = ene_cou_grad.ene
        state.epot = state.ene[state.instate]
        if self.coupling == "nacs":
            state.u = ene_cou_grad.u
            state.nac = ene_cou_grad.nac
            state.vk = self.vk_coupl_matrix(state.nac,state.vel)
            state.rho = propagator.elec_density(state)
        elif self.coupling == "wf_overlap":
            state.u = ene_cou_grad.u
            state.vk = ene_cou_grad.wf_ov 
            state.rho = propagator.elec_density(state)
        elif self.coupling == "non_coup":
            state.u = None 
            state.vk = None 
        state.ekin = self.cal_ekin(state.mass, state.vel)
        return grad_old

    def surface_hopping(self, state, ene_cou_grad, probs):            
        hop = None
        att = None
        succ = None
        aleatory = uniform(0,1)
        if self.coupling != "non_coup":
            acc_probs = cumsum(probs)
            total = sum(acc_probs)
            if total > 1.0:
                acc_probs /= total
            hopps = less(aleatory, acc_probs)
            if any(hopps):
                for i in range(self.nstates):
                    if hopps[i]:
                        state_new = state.states[i]
                        break
                #else:
                #    state_new = state.instate
                rescale = RescaleVelocity(state, ene_cou_grad)
                hop =  rescale.rescale_velocity(state, state_new) 
                if hop == "not":
                    state_new = state.instate
                    att = "yes"
                    succ = "not"
                else:
                    att = "yes"
                    succ = "yes"
            else:
                state_new = state.instate
            state.instate = state_new
            sur_hop = namedtuple("sur_hop", "aleatory acc_probs state_new hop att succ")
            return sur_hop(aleatory, acc_probs[state_new], state_new, hop, att, succ)
        elif self.coupling == "non_coup":
            iselected = probs.iselected 
            prob = probs.prob
            if prob > aleatory and iselected is not None:
                state_new = iselected
                rescale = RescaleVelocity(state, ene_cou_grad)
                hop =  rescale.rescale_velocity(state, iselected) 
                if hop == "not":
                    state_new = state.instate
            else:
                state_new = state.instate
            state.instate = state_new
            sur_hop = namedtuple("sur_hop", "aleatory prob state_new")
            return sur_hop(aleatory, prob, state_new)
            

    def new_surface(self, state, results, crd_new, t, dt):
        ene_cou_grad = self.get_ene_cou_grad(crd_new, state.instate)
        if self.coupling == "non_coup":
            vk_old = None 
            self.vel_old = state.vel
            propagator = Propagator(state)
            old_state = state.instate
            ise_prop = propagator.new_prob_grad(state, vk_old, ene_cou_grad, dt)
            #ise_prop = propagator.lz_select(state,dt)
            sur_hop = self.surface_hopping(state, ene_cou_grad, ise_prop.lz)
        elif self.coupling != "non_coup":
            vk_old = self.vk_coupl_matrix(state.nac,self.vel_old)
            self.vel_old = state.vel
            propagator = Propagator(state)
            grad_probs = propagator.new_prob_grad(state, vk_old, ene_cou_grad, dt)
            old_state = state.instate
            sur_hop = self.surface_hopping(state, ene_cou_grad, grad_probs.probs)
        state.ekin = self.cal_ekin(state.mass, state.vel)
        state.ene = ene_cou_grad.ene
        state.epot = state.ene[state.instate]
        results.save_db(t,state) #save variables in database
        results.print_var(t, dt, sur_hop, state) #printing variables 
        if self.coupling == "non_coup":
            ene_cou_grad = self.get_ene_cou_grad(crd_new, state.instate)
            return ene_cou_grad.grad        
        state.u = ene_cou_grad.u
        propagator.new_ncoeff(state, grad_probs, sur_hop.hop, sur_hop.att, sur_hop.succ)
        if self.coupling == "nacs":
            state.nac = ene_cou_grad.nac
            state.vk = self.vk_coupl_matrix(state.nac,state.vel)
        elif self.coupling == "wf_overlap":
            state.vk = ene_cou_grad.wf_ov 
        if old_state == sur_hop.state_new:
            return grad_probs.grad_new
        else:
            u_grad = self.get_hopp_u_grad(crd_new, state)
            propagator = Propagator(state)
            grad_new = propagator.new_hopp_grad(state, u_grad)
            return grad_new

class State(Colt):

    _user_input = """ 
    # chosen parameters
    db_file = :: existing_file 
    t = 0.0 :: float
    dt = 1.0 :: float
    mdsteps = 40000 :: float
    # instate is the initial state: 0 = G.S, 1 = E_1, ...
    instate = 1 :: int
    nstates = 2 :: int
    states = 0 1 :: ilist
    ncoeff = 0.0 1.0 :: flist
    # diagonal probability is not working yet
    prob = tully :: str :: tully, lz     
    rescale_vel = momentum :: str :: momentum, nacs 
    coupling = nacs :: str :: nacs, wf_overlap, non_coup
    method = Surface_Hopping :: str :: Surface_Hopping, Born_Oppenheimer  
    decoherence = EDC :: str :: EDC, IDC_A, IDC_S, No_DC 
    """
    
    def __init__(self, crd, vel, mass, model, t, dt, mdsteps, instate, nstates, states, ncoeff, prob, rescale_vel, coupling, method, decoherence, atomids):
        self.crd = crd
        self.natoms = len(crd)
        self.atomids = atomids
        self.vel = vel
        self.mass = mass
        if model == 1:
            self.model = True
        else:
            self.model = False
        self.t = t
        self.dt = dt
        self.mdsteps = mdsteps
        self.instate = instate
        self.nstates = nstates
        self.states = states
        self.ncoeff = ncoeff
        self.prob = prob
        self.rescale_vel = rescale_vel
        self.coupling = coupling
        if self.rescale_vel == "nacs" and self.coupling != "nacs":
            raise SystemExit("Wrong coupling method or wrong rescaling velocity approach")
        self.method = method
        self.decoherence = decoherence
        self.e_curr = None
        self.e_prev_step = None
        self.e_two_prev_steps = None
        self.ekin = 0
        self.epot = 0
        self.nac = {}
        self.ene = []
        self.vk = []
        self.u = []
        self.rho = []
        if isscalar(self.mass):
            self.natoms = 1
        elif isinstance(self.mass, ndarray) != True:
            self.natoms = array([self.mass])

    @classmethod
    def from_config(cls, config):
        crd, vel, mass, atomids, model = cls.read_db(config["db_file"])
        t = config['t']
        dt = config['dt']
        mdsteps = config['mdsteps']
        instate = config['instate']
        nstates = config['nstates']
        states = config['states']
        ncoeff = config['ncoeff']
        prob = config['prob']
        rescale_vel = config['rescale_vel']
        coupling = config['coupling']
        method = config['method']
        decoherence = config['decoherence']
        return cls(crd, vel, mass, model, t, dt, mdsteps, instate, nstates, states, ncoeff, prob, rescale_vel, coupling, method, decoherence, atomids)  

    @staticmethod
    def read_db(db_file):
        db = PySurfDB.load_database(db_file, read_only=True)
        crd = copy(db['crd'][0])
        vel = copy(db['veloc'][0])
        atomids = copy(db['atomids'])
        mass = copy(db['masses'])
        model = copy(db['model'])
        if model == 1:
            model = True
        else:
            model = False
        return crd, vel, mass, atomids, model

    @classmethod
    def from_initial(cls, crd, vel, mass, model, t, dt, mdsteps, instate, nstates, states, ncoeff, prob, rescale_vel, coupling, method, decoherence, atomids):
        return cls(crd, vel, mass, model, t, dt, mdsteps, instate, nstates, states, ncoeff, prob, rescale_vel, coupling, method, decoherence, atomids)

class PrintResults:
 
    def __init__(self):
        self.large = 110
        self.large_bo = 108
        self.dash = '-' * self.large
        self.dash_bo = '-' * self.large_bo
        self.gen_results = open("gen_results.out", "w")
        self.hopping = []
        self.tra_time = time()

    def norm_coeff(self, ncoeff):
        if isinstance(ncoeff, ndarray) != True:
            ncoeff = array(ncoeff)  
        return diag(outer(ncoeff, ncoeff.conj()).real)

    def save_db(self, t, state):
        nstates = state.nstates
        model = state.model
        nmodes = len(state.mass)
        prob = state.prob
        if isscalar(state.crd):
            natoms = 1
        else:
            natoms = len(state.crd)
        if state.method == "Surface_Hopping" and prob == "tully":
            if model:
                db = PySurfDB.generate_database("results.db", data=["crd","veloc","energy","time","ekin","epot","etot","fosc","currstate"], dimensions ={"nmodes":nmodes, "nstates":nstates}, model = model)
            else:
                db = PySurfDB.generate_database("results.db", data=["crd","veloc","energy","time","ekin","epot","etot","fosc","currstate"], dimensions ={"natoms":natoms, "nstates":nstates}, model = model)
            db.set("currstate",state.instate)
            db.set("fosc",self.norm_coeff(state.ncoeff))
        elif state.method == "Surface_Hopping" and prob == "lz":
            if model:
                db = PySurfDB.generate_database("results.db", data=["crd","veloc","energy","time","ekin","epot","etot","currstate"], dimensions ={"nmodes":nmodes, "nstates":nstates}, model = model)
            else:
                db = PySurfDB.generate_database("results.db", data=["crd","veloc","energy","time","ekin","epot","etot","currstate"], dimensions ={"natoms":natoms, "nstates":nstates}, model = model)
            db.set("currstate",state.instate)
        elif state.method == "Born_Oppenheimer":
            db = PySurfDB.generate_database("results.db", data=["crd","veloc","energy","time","ekin","epot","etot"], dimensions ={"natoms":natoms, "nstates":nstates}, model = model)
        db.set("crd",state.crd)
        db.set("veloc",state.vel)
        db.set("energy",state.ene)
        db.set("time",t)
        db.set("ekin",state.ekin)
        db.set("epot",state.epot)
        db.set("etot",state.ekin+state.epot)
        db.increase  # It increases the frame

    def dis_dimer(self, a,b):
        return sqrt(sum((a-b)**2))
    
    def print_acknowledgment(self, state):
        title = " Trajectory Surface Hopping Module "
        based = " This module uses the tools implemented in PySurf"
        contributors = " Module implemented by: Edison Salazar, Maximilian Menger, and Shirin Faraji "
        vel = state.vel
        crd = state.crd
        prob = state.prob
        coupling = state.coupling
        rescale_vel = state.rescale_vel 
        dt = state.dt
        mdsteps = state.mdsteps
        instate = state.instate
        self.instate = instate
        nstates = state.nstates
        ncoeff = self.norm_coeff(state.ncoeff)
        decoherence = state.decoherence
        ack = namedtuple("ack", "title vel crd based actors prob coupling rescale_vel dt mdsteps instate nstates ncoeff decoherence")
        return ack(title, vel, crd, based, contributors, prob, coupling, rescale_vel, dt, mdsteps, instate, nstates, ncoeff, decoherence)

    def print_head(self, state):
        ack = self.print_acknowledgment(state)  
        if state.method == "Surface_Hopping":
            self.gen_results.write(f"\n{ack.title:=^{self.large}}\n")
            self.gen_results.write(f"\n{ack.based:^{self.large}}\n")
            self.gen_results.write(f"{ack.actors:^{self.large}}\n")        
            self.gen_results.write(f"\nInitial parameters:\n")
            self.gen_results.write(f"   Time step: {ack.dt}\n")
            self.gen_results.write(f"   MD steps: {ack.mdsteps}\n")
            self.gen_results.write(f"   Number of states: {ack.nstates}\n")
            self.gen_results.write(f"   Initial population: {ack.ncoeff}\n")
            self.gen_results.write(f"   Initial state: {ack.instate}\n")
            self.gen_results.write(f"   Probability method: {ack.prob}\n")
            self.gen_results.write(f"   Coupling: {ack.coupling}\n")
            self.gen_results.write(f"   Rescale of velocity: {ack.rescale_vel}\n")
            self.gen_results.write(f"   Decoherence: {ack.decoherence}\n")
            self.gen_results.write(f"Computing a trajectory surface hopping simulation:\n")
            self.gen_results.write(self.dash + "\n")
            head = namedtuple("head","steps t ekin epot etotal hopp random state")
            head = head("MD_steps", "Time", "E_kinetic", "E_potential", "E_total", "Hopping_P", "Random", "State")
            self.gen_results.write(f"{head.steps:>10s} {head.t:>10s} {head.ekin:>15s} {head.epot:>17s} {head.etotal:>13s}"\
                    f"{head.hopp:>15s} {head.random:>11s} {head.state:>11s} \n")
            self.gen_results.write(self.dash + "\n")
        elif state.method == "Born_Oppenheimer":
            self.t_crd_vel_ene_popu = open("t_crd_vel_ene_popu.csv", "w")
            self.gen_results.write(f"\n{ack.title:=^{self.large_bo}}\n")
            self.gen_results.write(f"\n{ack.based:^{self.large_bo}}\n")
            self.gen_results.write(f"{ack.actors:^{self.large_bo}}\n")        
            self.gen_results.write(f"Initial parameters:\n")
            self.gen_results.write(f"   Initial position:\n")
            self.gen_results.write(f"   {self.dis_dimer(ack.crd[0],ack.crd[1]):>0.4f}\n")
            self.gen_results.write(f"   Initial velocity:\n")
            self.gen_results.write(f"   {self.dis_dimer(ack.vel[0],ack.vel[1]):>0.4f}\n")
            self.gen_results.write(f"   Time step: {ack.dt}\n")
            self.gen_results.write(f"   MD steps: {ack.mdsteps}\n")
            self.gen_results.write(f"   Active state: {ack.instate}\n")
            self.gen_results.write(f"Computing a Born Oppenheimer simulation:\n")
            self.gen_results.write(self.dash_bo + "\n")
            head = namedtuple("head","steps t dis dis_vel ekin epot etotal state")
            head = head("MD_steps", "Time", "D_r1-r2", "D_v1-v2", "E_kinetic", "E_potential", "E_total", "State")
            self.gen_results.write(f"{head.steps:>10s} {head.t:>10s} {head.dis:>12s} {head.dis_vel:>12s}"\
                    f"{head.ekin:>15s} {head.epot:>17s} {head.etotal:>13s} {head.state:>11s} \n")
            self.gen_results.write(self.dash_bo + "\n")
            self.t_crd_vel_ene_popu.write(f"{head.t},{head.dis},{head.dis_vel},{head.ekin},"\
                    f"{head.epot},{head.etotal},{head.state}\n")

    def print_var(self, t, dt, sur_hop, state):        
        var = namedtuple("var","steps t ekin epot etotal hopp random state")
        if state.prob == "tully":
            var = var(int(t/dt),t,state.ekin,state.epot,state.ekin + state.epot,\
                    sur_hop.acc_probs,sur_hop.aleatory,state.instate)
        elif state.prob == "lz":
            var = var(int(t/dt),t,state.ekin,state.epot,state.ekin + state.epot,\
                    sur_hop.prob,sur_hop.aleatory,state.instate)
        self.gen_results.write(f"{var.steps:>8.0f} {var.t:>12.2f} {var.ekin:>15.3f} {var.epot:>17.4f}"\
                    f"{var.etotal:>13.4f} {var.hopp:>15.5f} {var.random:>11.5f} {var.state:>11.0f}\n")
        if var.state != self.instate:
            self.hopping.append(f"Hopping from state {self.instate} to state {state.instate}"\
                                f" in step: {var.steps}, at the time step: {var.t}")
            self.instate = var.state

    def print_bh_var(self, t, dt, state):        
        var = namedtuple("var","steps t dis dis_vel ekin epot etotal state")
        var = var(int(t/dt),t,self.dis_dimer(state.crd[0],state.crd[1]),self.dis_dimer(state.vel[0],\
                    state.vel[1]),state.ekin,state.epot,state.ekin + state.epot,state.instate)
        self.gen_results.write(f"{var.steps:>8.0f} {var.t:>12.2f} {var.dis:>12.4f} {var.dis_vel:>12.4f}"\
                    f"{var.ekin:>15.3f} {var.epot:>17.4f} {var.etotal:>13.4f} {var.state:>11.0f}\n")
        self.t_crd_vel_ene_popu.write(f"{var.t:>0.3f},{var.dis:>0.8f},{var.dis_vel:>0.8f},"\
                    f"{var.ekin:>0.8f},{var.epot:>0.8f},{var.etotal:>0.8f},{var.state:>0.0f}\n")

    def print_bottom(self, state):
        if state.method == "Surface_Hopping":
            self.gen_results.write(self.dash + "\n")
            self.gen_results.write(f"Some important variables are printed in results.db\n")
            if self.hopping:
                for i in range(len(self.hopping)):
                    self.gen_results.write(f"{self.hopping[i]}\n")
                #if i == 0:
                #    print(f"There is {i+1} hop")
                #else:
                #    print(f"There are {i+1} hops")
            else:
                self.gen_results.write(f"No hoppings achieved\n")
        elif state.method == "Born_Oppenheimer": 
            self.t_crd_vel_ene_popu.close()
            self.gen_results.write(self.dash_bo + "\n")
            self.gen_results.write(f"Some important variables are printed in t_crd_vel_ene_popu.csv and results.db\n")
        time_seg = time()-self.tra_time
        day = time_seg // (24*3600)
        time_seg = time_seg % (24*3600)
        hour = time_seg // 3600
        time_seg %= 3600
        minutes = time_seg // 60
        time_seg %= 60
        seconds = time_seg
        self.gen_results.write(f"Total job time: {day:>0.0f}:{hour:>0.0f}:{minutes:>0.0f}:{seconds:>0.0f}\n")
        
        self.gen_results.write(f"{ctime()}")
        self.gen_results.close()

if __name__=="__main__":
    elec_state = State.from_questions(config = "prop.inp")
    DY = VelocityVerletPropagator(elec_state)    
    try:
        result_2 = DY.run()
    except SystemExit as err:
        print("An error:", err) 

    
