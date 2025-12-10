package NitrogenMembranePDE  /**
         * PDE-based membrane model (1D diffusion across thickness) for 2 components:
         * O2 and N2.
         */

  model MembranePDE2Comp
    // ----------- Connectors -----------
    Modelica.Blocks.Interfaces.RealInput p_feed "Feed pressure [Pa]";
    Modelica.Blocks.Interfaces.RealOutput yN2_ret "N2 mole fraction in retentate";
    Modelica.Blocks.Interfaces.RealOutput yO2_ret "O2 mole fraction in retentate";
    Modelica.Blocks.Interfaces.RealOutput F_ret "Retentate molar flow [mol/s]";
    Modelica.Blocks.Interfaces.RealOutput Fperm_O2 "O2 permeate flow [mol/s]";
    Modelica.Blocks.Interfaces.RealOutput Fperm_N2 "N2 permeate flow [mol/s]";
  
    // ----------- Feed conditions -----------
    parameter Real F_in   = 1.0  "Feed molar flow [mol/s]";
    parameter Real yO2_in = 0.21 "O2 mole fraction in feed";
    parameter Real yN2_in = 0.79 "N2 mole fraction in feed";
  
    // ----------- Permeate side conditions -----------
    parameter Real p_perm   = 1e5 "Permeate pressure [Pa]";
    parameter Real yO2_perm = 0.0 "O2 mole fraction on permeate side";
    parameter Real yN2_perm = 0.0 "N2 mole fraction on permeate side";
  
    // ----------- Geometry -----------
    parameter Real A     = 1.0   "Membrane area [m2]";      // << smaller, more realistic
    parameter Real delta = 2e-7  "Membrane thickness [m]";
  
    // ----------- Transport properties -----------
    parameter Real D_O2 = 4e-10   "O2 diffusivity [m2/s]";
    parameter Real D_N2 = 0.6e-10 "N2 diffusivity [m2/s]";
    parameter Real S_O2 = 2e-3    "O2 solubility [mol/(m3*Pa)]";
    parameter Real S_N2 = 0.6e-3  "N2 solubility [mol/(m3*Pa)]";
  
    // ----------- Discretisation -----------
    parameter Integer N(min = 3) = 10 "Number of spatial nodes";
    parameter Real dx = delta/(N - 1) "Spatial step [m]";
  
    // Concentration profiles in the polymer
    Real C_O2[N] "O2 concentration inside membrane [mol/m3]";
    Real C_N2[N] "N2 concentration inside membrane [mol/m3]";
  
    // Boundary partial pressures
    Real pO2_feed, pN2_feed;
    Real pO2_perm, pN2_perm;
  
    // Fluxes [mol/(m2*s)]
    Real J_O2, J_N2;
  
    // ---- NEW: component flows ----
    Real Fperm_O2_raw, Fperm_N2_raw;
    Real Fperm_O2_lim, Fperm_N2_lim;
    Real F_O2_ret, F_N2_ret;
  equation
    // ---- Partial pressures on feed and permeate sides ----
    pO2_feed = yO2_in*p_feed;
    pN2_feed = yN2_in*p_feed;
    pO2_perm = yO2_perm*p_perm;
    pN2_perm = yN2_perm*p_perm;
  
    // ---- Boundary conditions (Henry's law) ----
    C_O2[1] = S_O2*pO2_feed;
    C_N2[1] = S_N2*pN2_feed;
    C_O2[N] = S_O2*pO2_perm;
    C_N2[N] = S_N2*pN2_perm;
  
    // ---- Internal nodes: diffusion PDE -> ODEs ----
    for j in 2:(N - 1) loop
      der(C_O2[j]) = D_O2*(C_O2[j + 1] - 2*C_O2[j] + C_O2[j - 1])/dx^2;
      der(C_N2[j]) = D_N2*(C_N2[j + 1] - 2*C_N2[j] + C_N2[j - 1])/dx^2;
    end for;
  
    // ---- Fluxes at permeate side ----
    J_O2 = -D_O2*(C_O2[N] - C_O2[N - 1])/dx;
    J_N2 = -D_N2*(C_N2[N] - C_N2[N - 1])/dx;
  
    // ---- Raw permeate molar flows from PDE flux ----
    Fperm_O2_raw = A*J_O2;
    Fperm_N2_raw = A*J_N2;
  
    // ---- Limit permeation: 0 <= Fperm_i <= 0.99 * feed component flow ----
    Fperm_O2_lim = min(max(Fperm_O2_raw, 0), 0.99*F_in*yO2_in);
    Fperm_N2_lim = min(max(Fperm_N2_raw, 0), 0.99*F_in*yN2_in);
  
    // Exposed permeate flows
    Fperm_O2 = Fperm_O2_lim;
    Fperm_N2 = Fperm_N2_lim;
  
    // ---- Retentate component flows ----
    F_O2_ret = F_in*yO2_in - Fperm_O2_lim;
    F_N2_ret = F_in*yN2_in - Fperm_N2_lim;
  
    // ---- Total retentate flow and composition ----
    F_ret   = F_O2_ret + F_N2_ret;
    yO2_ret = F_O2_ret/F_ret;
    yN2_ret = F_N2_ret/F_ret;
  end MembranePDE2Comp;

  model ReceiverTank_Block
    // --- Connectors ---
    Modelica.Blocks.Interfaces.RealInput  F_in  "N2 inflow [mol/s]";
    Modelica.Blocks.Interfaces.RealInput  F_out "N2 outflow [mol/s]";
    Modelica.Blocks.Interfaces.RealOutput p     "Tank pressure [Pa]";
    Modelica.Blocks.Interfaces.RealOutput n     "Moles of N2 in tank [mol]";
  
    // --- Parameters ---
    parameter Real V = 0.05 "Tank volume [m3]";
    parameter Real T = 293  "Tank temperature [K]";
    constant  Real R = 8.314 "Gas constant [J/(mol*K)]";
  
    // --- Net flow and integrator ---
    Modelica.Blocks.Math.Add add(k1=+1, k2=-1) 
      "F_in - F_out";
    Modelica.Blocks.Continuous.Integrator int(
      k       = 1,
      y_start = 0) 
      "Integrates net flow to get moles";
  equation
    // Net flow
    connect(F_in,  add.u1);
    connect(F_out, add.u2);
  
    // Integrate net flow
    connect(add.y, int.u);
  
    // Outputs
    n = int.y;
    p = n*R*T / V;
  end ReceiverTank_Block;

  model SimpleValve
    Modelica.Blocks.Interfaces.RealInput  p_tank "Tank pressure [Pa]";
    Modelica.Blocks.Interfaces.RealOutput F_out  "N2 outflow [mol/s]";
  
    parameter Real p_atm = 1e5  "Atmospheric pressure [Pa]";
    parameter Real Kv    = 1e-7 "Valve coefficient [mol/(s*Pa)]";
  equation
    // Linear valve: outflow proportional to overpressure above atm
    F_out = if p_tank > p_atm then Kv*(p_tank - p_atm) else 0;
  end SimpleValve;

  model NitrogenSystemPDE
    import Modelica.Blocks.Sources;
  
    // Pump: feed pressure step
    Sources.Step pumpP(
      startTime = 5,
      offset    = 1e5,   // 5 bar
      height    = 0.5e5)   // +3 bar -> 8 bar
      "Represents compressor inlet pressure";
  
    // PDE-based membrane
    MembranePDE2Comp mem;
  
    // Receiver tank (ideal gas + integration)
    ReceiverTank_Block tank(
      V = 0.05,
      T = 293);
  
    // Valve controlled by tank pressure
    SimpleValve valve(
      Kv = 1e-7);
  equation
    // Pump → membrane
    connect(pumpP.y, mem.p_feed);
  
    // Membrane retentate (mostly N2) → tank inflow
    tank.F_in = mem.F_ret * mem.yN2_ret;
  
    // Tank pressure → valve, valve outflow → tank
    connect(tank.p,      valve.p_tank);
    connect(valve.F_out, tank.F_out);
  end NitrogenSystemPDE;
end NitrogenMembranePDE;
