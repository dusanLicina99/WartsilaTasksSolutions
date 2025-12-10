model VentilationRoom
  // Use built-in fluid & thermal components
  inner Modelica.Fluid.System system;

  // -------- PARAMETERS --------
  parameter Real V         = 50      "Room volume [m3]";
  parameter Real Qh        = 1500    "Heater power [W]";
  parameter Real mIn       = 0.05    "Ventilation inflow [kg/s]";
  parameter Real mLeak     = 0.002   "Base leak (door closed) [kg/s]";
  parameter Real mDoorMax  = 0.15    "Extra leak when door is open [kg/s]";
  parameter Real tDoorOpen = 10      "Door opens at this time [s]";

  parameter Real T_supply  = 298.15  "Ventilation supply temperature (25°C)";
  parameter Real T_out     = 273.15  "Outside temperature (0°C)";

  // -------- STATE: gas / old-air fraction --------
  Real Yg(start = 1)       "Old air / contaminant mass fraction";
  Real mLeakEff            "Effective leak used in Yg equation [kg/s]";

  // -------- FLUID COMPONENTS --------
  // Room volume (3 ports: vent in, leak in, outlet)
  Modelica.Fluid.Vessels.ClosedVolume room(
    V               = V,
    nPorts          = 3,
    use_portsData   = false,
    redeclare package Medium = Modelica.Media.Air.SimpleAir);

  // Heater
  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow heater(
    Q_flow = Qh,
    alpha  = 500);

  // Ventilation inflow (fresh air at T_supply)
  Modelica.Fluid.Sources.MassFlowSource_T ventIn(
    m_flow        = mIn,
    T             = T_supply,
    nPorts        = 1,
    redeclare package Medium = Modelica.Media.Air.SimpleAir);

  // Leak / door inflow: mass flow controlled via m_flow_in
  Modelica.Fluid.Sources.MassFlowSource_T leakIn(
    use_m_flow_in = true,
    m_flow        = mLeak,        // nominal value (not used when use_m_flow_in = true)
    T             = T_out,
    nPorts        = 1,
    redeclare package Medium = Modelica.Media.Air.SimpleAir);

  // Outlet to outside (pressure boundary)
  Modelica.Fluid.Sources.Boundary_pT outlet(
    p       = 101325,
    T       = T_out,
    nPorts  = 1,
    redeclare package Medium = Modelica.Media.Air.SimpleAir);

  // Step signal: at tDoorOpen leak increases by mDoorMax
  Modelica.Blocks.Sources.Step leakFlow(
    height    = mDoorMax,
    offset    = mLeak,
    startTime = tDoorOpen);

equation
  // Fluid connections
  connect(heater.port,      room.heatPort);
  connect(ventIn.ports[1],  room.ports[1]);
  connect(leakIn.ports[1],  room.ports[2]);
  connect(outlet.ports[1],  room.ports[3]);

  // Control leak mass-flow with step block (scalar -> scalar)
  connect(leakFlow.y, leakIn.m_flow_in);

  // Effective leak used in our simple Yg model
  mLeakEff = leakFlow.y;

  // Gas composition ODE (rho ≈ 1.2 kg/m3, well-mixed)
  der(Yg) = (mLeakEff - Yg*(mIn + mLeakEff)) / (1.2*V);

end VentilationRoom;


















