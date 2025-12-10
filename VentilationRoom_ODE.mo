model VentilationRoom_ODE
  // -------- PARAMETERS --------
  parameter Real V         = 50      "Room volume [m3]";
  parameter Real rho       = 1.2     "Air density [kg/m3]";
  parameter Real cp        = 1000    "Specific heat of air [J/(kg*K)]";

  parameter Real Qh        = 1500    "Heater power [W]";
  parameter Real mIn       = 0.05    "Ventilation inflow [kg/s]";
  parameter Real mLeak     = 0.002   "Base leak outflow [kg/s]";
  parameter Real mDoorMax  = 0.20    "Door leak when door is open [kg/s]";

  parameter Real T_out     = 273.15  "Outside air temperature [K] (0°C)";
  parameter Real T_supply  = 298.15  "Ventilation supply temperature [K] (25°C)";

  parameter Real tDoorOpen = 10      "Time when door opens [s]";

  // -------- DERIVED --------
  Real m "Mass of air in room [kg]";
  Real doorOpen "0 = closed, 1 = open";
  Real mLeakTotal "Total leak outflow [kg/s]";

  // -------- STATES --------
  Real T(start = 293)   "Room temperature [K]";
  Real Yg(start = 1)    "Old air / contaminant mass fraction [-]";

  // Optional: for plotting in ppm
  Real Cppm "Old air concentration [ppm]";

equation
  // Total mass of air (assumed constant)
  m = rho * V;

  // Door logic: closed before tDoorOpen, open after
  doorOpen = if time < tDoorOpen then 0 else 1;

  // Total leak outflow = base leak + door leak
  mLeakTotal = mLeak + doorOpen * mDoorMax;

  // ---------- ENERGY BALANCE (temperature dynamics) ----------
  // Qh: heater power [W]
  // mIn*cp*(T_supply - T): ventilation at T_supply
  // mLeakTotal*cp*(T_out - T): leaks replaced by outside air at T_out
  der(T) =
    ( Qh
    + mIn*cp*(T_supply - T)
    + mLeakTotal*cp*(T_out    - T)
    ) / (m*cp);

  // ---------- OLD-AIR / GAS FRACTION DYNAMICS ----------
  // Fresh air (mIn) + leaks (mLeakTotal) flush out old air
  der(Yg) = - Yg * (mIn + mLeakTotal) / m;

  // ppm for nicer plots
  Cppm = Yg * 1e6;

end VentilationRoom_ODE;
