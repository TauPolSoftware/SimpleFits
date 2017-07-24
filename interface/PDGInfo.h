#ifndef PDGInfo_h
#define PDGInfo_h

#include <sstream>

class PDGInfo {
 public:
  static double pi_mass(){return  0.13957018;}
  static double tau_mass(){return 1.77682;}
  static double nu_mass(){return  0.0;}

  static double pi_mass_MCGen(){return  0.139;}
  static double tau_mass_MCGen(){return 1.777;}
  static double nu_mass_MCGen(){return  0.0;}

  enum PDGMCNumbering {
    d = 1 ,
    anti_d = -1 ,
    u = 2 ,
    anti_u = -2 ,
    s = 3 ,
    anti_s = -3 ,
    c = 4 ,
    anti_c = -4 ,
    b = 5 ,
    anti_b = -5 ,
    t = 6 ,
    anti_t = -6 ,
    l = 7 ,
    anti_l = -7 ,
    h = 8 ,
    anti_h = -8 ,
    g = 21 ,
    e_minus = 11 ,
    e_plus = -11 ,
    nu_e = 12 ,
    anti_nu_e = -12 ,
    mu_minus = 13 ,
    mu_plus = -13 ,
    nu_mu = 14 ,
    anti_nu_mu = -14 ,
    tau_minus = 15 ,
    tau_plus = -15 ,
    nu_tau = 16 ,
    anti_nu_tau = -16 ,
    L_minus = 17 ,
    L_plus = -17 ,
    nu_L = 18 ,
    anti_nu_L = -18 ,
    gamma = 22 ,
    Z0 = 23 ,
    W_plus = 24 ,
    W_minus = -24 ,
    Higgs0 = 25 ,
    reggeon = 28 ,
    pomeron = 29 ,
    Z_prime0 = 32 ,
    Z_prime_prime0 = 33 ,
    W_prime_plus = 34 ,
    W_prime_minus = -34 ,
    Higgs_prime0 = 35 ,
    A0 = 36 ,
    Higgs_plus = 37 ,
    Higgs_minus = -37 ,
    pi0 = 111 ,
    pi_plus = 211 ,
    pi_minus = -211 ,
    pi_diffr_plus = 210 ,
    pi_diffr_minus = -210 ,
    pi_2S0 = 20111 ,
    pi_2S_plus = 20211 ,
    pi_2S_minus = -20211 ,
    eta = 221 ,
    eta_2S = 20221 ,
    eta_prime = 331 ,
    rho0 = 113 ,
    rho_plus = 213 ,
    rho_minus = -213 ,
    rho_2S0 = 30113 ,
    rho_2S_plus = 30213 ,
    rho_2S_minus = -30213 ,
    rho_3S0 = 40113 ,
    rho_3S_plus = 40213 ,
    rho_3S_minus = -40213 ,
    omega = 223 ,
    omega_2S = 30223 ,
    phi = 333 ,
    a_00 = 10111 ,
    a_0_plus = 10211 ,
    a_0_minus = -10211 ,
    f_0 = 10221 ,
    f_prime_0 = 10331 ,
    b_10 = 10113 ,
    b_1_plus = 10213 ,
    b_1_minus = -10213 ,
    h_1 = 10223 ,
    h_prime_1 = 10333 ,
    a_10 = 20113 ,
    a_1_plus = 20213 ,
    a_1_minus = -20213 ,
    f_1 = 20223 ,
    f_prime_1 = 20333 ,
    a_20 = 115 ,
    a_2_plus = 215 ,
    a_2_minus = -215 ,
    f_2 = 225 ,
    f_prime_2 = 335 ,
    K0 = 311 ,
    anti_K0 = -311 ,
    K_S0 = 310 ,
    K_L0 = 130 ,
    K_plus = 321 ,
    K_minus = -321 ,
    K_star0 = 313 ,
    anti_K_star0 = -313 ,
    K_star_plus = 323 ,
    K_star_minus = -323 ,
    K_0_star0 = 10311 ,
    anti_K_0_star0 = -10311 ,
    K_0_star_plus = 10321 ,
    K_0_star_minus = -10321 ,
    K_10 = 10313 ,
    anti_K_10 = -10313 ,
    K_1_plus = 10323 ,
    K_1_minus = -10323 ,
    K_2_star0 = 315 ,
    anti_K_2_star0 = -315 ,
    K_2_star_plus = 325 ,
    K_2_star_minus = -325 ,
    K_prime_10 = 20313 ,
    anti_K_prime_10 = -20313 ,
    K_prime_1_plus = 20323 ,
    K_prime_1_minus = -20323 
   };

  static std::string pdgIdToName(int pdgId) {
	   if(pdgId == 1)        return "d";
	   if(pdgId == -1)       return "anti_d";
	   if(pdgId == 2)        return "u";
	   if(pdgId == -2)       return "anti_u";
	   if(pdgId == 3)        return "s";
	   if(pdgId == -3)       return "anti_s";
	   if(pdgId == 4)        return "c";
	   if(pdgId == -4)       return "anti_c";
	   if(pdgId == 5)        return "b";
	   if(pdgId == -5)       return "anti_b";
	   if(pdgId == 6)        return "t";
	   if(pdgId == -6)       return "anti_t";
	   if(pdgId == 7)        return "l";
	   if(pdgId == -7)       return "anti_l";
	   if(pdgId == 8)        return "h";
	   if(pdgId == -8)       return "anti_h";
	   if(pdgId == 21)       return "g";
	   if(pdgId == 11)       return "e-";
	   if(pdgId == -11)      return "e+";
	   if(pdgId == 12)       return "nu_e";
	   if(pdgId == -12)      return "anti_nu_e";
	   if(pdgId == 13)       return "mu-";
	   if(pdgId == -13)      return "mu+";
	   if(pdgId == 14)       return "nu_mu";
	   if(pdgId == -14)      return "anti_nu_mu";
	   if(pdgId == 15)       return "tau-";
	   if(pdgId == -15)      return "tau+";
	   if(pdgId == 16)       return "nu_tau";
	   if(pdgId == -16)      return "anti_nu_tau";
	   if(pdgId == 17)       return "L-";
	   if(pdgId == -17)      return "L+";
	   if(pdgId == 18)       return "nu_L";
	   if(pdgId == -18)      return "anti_nu_L";
	   if(pdgId == 22)       return "gamma";
	   if(pdgId == 23)       return "Z0";
	   if(pdgId == 24)       return "W+";
	   if(pdgId == -24)      return "W-";
	   if(pdgId == 25)       return "Higgs0";
	   if(pdgId == 28)       return "reggeon";
	   if(pdgId == 29)       return "pomeron";
	   if(pdgId == 32)       return "Z_prime0";
	   if(pdgId == 33)       return "Z_prime_prime0";
	   if(pdgId == 34)       return "W_prime+";
	   if(pdgId == -34)      return "W_prime-";
	   if(pdgId == 35)       return "Higgs_prime0";
	   if(pdgId == 36)       return "A0";
	   if(pdgId == 37)       return "Higgs+";
	   if(pdgId == -37)      return "Higgs-";
	   if(pdgId == 111)      return "pi0";
	   if(pdgId == 211)      return "pi+";
	   if(pdgId == -211)     return "pi-";
	   if(pdgId == 210)      return "pi_diffr+";
	   if(pdgId == -210)     return "pi_diffr-";
	   if(pdgId == 20111)    return "pi_2S0";
	   if(pdgId == 20211)    return "pi_2S+";
	   if(pdgId == -20211)   return "pi_2S-";
	   if(pdgId == 221)      return "eta";
	   if(pdgId == 20221)    return "eta_2S";
	   if(pdgId == 331)      return "eta_prime";
	   if(pdgId == 113)      return "rho0";
	   if(pdgId == 213)      return "rho+";
	   if(pdgId == -213)     return "rho-";
	   if(pdgId == 30113)    return "rho_2S0";
	   if(pdgId == 30213)    return "rho_2S+";
	   if(pdgId == -30213)   return "rho_2S-";
	   if(pdgId == 40113)    return "rho_3S0";
	   if(pdgId == 40213)    return "rho_3S+";
	   if(pdgId == -40213)   return "rho_3S-";
	   if(pdgId == 223)      return "omega";
	   if(pdgId == 30223)    return "omega_2S";
	   if(pdgId == 333)      return "phi";
	   if(pdgId == 10111)    return "a_00";
	   if(pdgId == 10211)    return "a_0+";
	   if(pdgId == -10211)   return "a_0-";
	   if(pdgId == 10221)    return "f_0";
	   if(pdgId == 10331)    return "f_prime_0";
	   if(pdgId == 10113)    return "b_10";
	   if(pdgId == 10213)    return "b_1+";
	   if(pdgId == -10213)   return "b_1-";
	   if(pdgId == 10223)    return "h_1";
	   if(pdgId == 10333)    return "h_prime_1";
	   if(pdgId == 20113)    return "a_10";
	   if(pdgId == 20213)    return "a_1+";
	   if(pdgId == -20213)   return "a_1-";
	   if(pdgId == 20223)    return "f_1";
	   if(pdgId == 20333)    return "f_prime_1";
	   if(pdgId == 115)      return "a_20";
	   if(pdgId == 215)      return "a_2+";
	   if(pdgId == -215)     return "a_2-";
	   if(pdgId == 225)      return "f_2";
	   if(pdgId == 335)      return "f_prime_2";
	   if(pdgId == 311)      return "K0";
	   if(pdgId == -311)     return "anti_K0";
	   if(pdgId == 310)      return "K_S0";
	   if(pdgId == 130)      return "K_L0";
	   if(pdgId == 321)      return "K+";
	   if(pdgId == -321)     return "K-";
	   if(pdgId == 313)      return "K_star0";
	   if(pdgId == -313)     return "anti_K_star0";
	   if(pdgId == 323)      return "K_star+";
	   if(pdgId == -323)     return "K_star-";
	   if(pdgId == 10311)    return "K_0_star0";
	   if(pdgId == -10311)   return "anti_K_0_star0";
	   if(pdgId == 10321)    return "K_0_star+";
	   if(pdgId == -10321)   return "K_0_star-";
	   if(pdgId == 10313)    return "K_10";
	   if(pdgId == -10313)   return "anti_K_10";
	   if(pdgId == 10323)    return "K_1+";
	   if(pdgId == -10323)   return "K_1-";
	   if(pdgId == 315)      return "K_2_star0";
	   if(pdgId == -315)     return "anti_K_2_star0";
	   if(pdgId == 325)      return "K_2_star+";
	   if(pdgId == -325)     return "K_2_star-";
	   if(pdgId == 20313)    return "K_prime_10";
	   if(pdgId == -20313)   return "anti_K_prime_10";
	   if(pdgId == 20323)    return "K_prime_1+";
	   if(pdgId == -20323)   return "K_prime_1-";

	   std::stringstream out;
	   out << "unknown ID = " << pdgId;
	   return out.str();
   };
};
#endif
