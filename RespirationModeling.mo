within ;
package RespirationModeling
  package Parts
    model Athmosphere "gaseous composition of the atmosphere"
      parameter Modelica.Units.SI.VolumeFraction FO2d = 0.2095; //O2 volume fraction in STPD
      parameter Modelica.Units.SI.VolumeFraction FCO2d = 0.0004;//CO2 volume fraction in STPD
      parameter Modelica.Units.SI.VolumeFraction FN2d = 0.7808;//N2 volume fraction in STPD
      parameter Modelica.Units.SI.VolumeFraction FArd = 0.0093;//Argon volume fraction in STPD
      parameter Modelica.Units.SI.MolarVolume VmO2stand=22.3922/1000 "O2 molar volume at 0°C 760 mmHg";
      parameter Modelica.Units.SI.MolarVolume VmCO2stand=22.2629/1000 "CO2 molar volume at 0°C 760 mmHg";
      parameter Modelica.Units.SI.MolarVolume VmN2stand= 22.4037/1000 "N2 molar volume at 0°C 760 mmHg";
      parameter Modelica.Units.SI.MolarVolume VmArstand=22.3929/1000 "Argon molar volume at 0°C 760 mmHg";
      parameter Modelica.Units.SI.Pressure Pstand= 101325; //pressure 760 mmHg in Pa
      Bodylight.Types.RealIO.PressureInput PB "pressure" annotation (Placement(
            transformation(extent={{-270,6},{-230,46}}), iconTransformation(extent={{-122,62},
                {-100,84}})));
      Bodylight.Types.RealIO.TemperatureInput temp "temperature" annotation (
          Placement(transformation(extent={{-270,2},{-230,42}}), iconTransformation(
              extent={{-122,10},{-100,32}})));
      Bodylight.Types.RealIO.FractionInput hum "relative air humidity"
        annotation (Placement(transformation(extent={{-272,-6},{-232,36}}),
            iconTransformation(extent={{-124,-46},{-100,-22}})));
      Bodylight.Types.RealIO.FractionOutput FO2 "O2 volume fraction"
        annotation (Placement(transformation(extent={{-260,-24},{-240,-4}}),
            iconTransformation(extent={{100,72},{120,92}})));
      Bodylight.Types.RealIO.ConcentrationOutput O2conc "O2 voncentration"
        annotation (Placement(transformation(extent={{-260,-8},{-240,12}}),
            iconTransformation(extent={{100,52},{120,72}})));
      Bodylight.Types.RealIO.FractionOutput FCO2 "CO2 volume fraction"
        annotation (Placement(transformation(extent={{-260,-24},{-240,-4}}),
            iconTransformation(extent={{100,12},{120,32}})));
      Bodylight.Types.RealIO.ConcentrationOutput CO2conc "CO2 voncentration"
        annotation (Placement(transformation(extent={{-260,-8},{-240,12}}),
            iconTransformation(extent={{100,-8},{120,12}})));
      Bodylight.Types.RealIO.FractionOutput FN2 "N2 volume fraction"
        annotation (Placement(transformation(extent={{-260,-24},{-240,-4}}),
            iconTransformation(extent={{100,-46},{120,-26}})));
      Bodylight.Types.RealIO.ConcentrationOutput N2conc "N2 voncentration"
        annotation (Placement(transformation(extent={{-260,-8},{-240,12}}),
            iconTransformation(extent={{100,-66},{120,-46}})));
      Bodylight.Types.RealIO.FractionOutput FAr "Ar volume fraction"
        annotation (Placement(transformation(extent={{-260,-24},{-240,-4}}),
            iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={44,-110})));
      Bodylight.Types.RealIO.ConcentrationOutput Arconc "Ar voncentration"
        annotation (Placement(transformation(extent={{-260,-8},{-240,12}}),
            iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={24,-110})));
      Bodylight.Types.RealIO.FractionOutput FH2O "vapor fraction" annotation (
          Placement(transformation(extent={{-260,-24},{-240,-4}}),
            iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-38,-110})));
      Bodylight.Types.RealIO.PressureOutput pO2 "O2 partial pressure "
        annotation (Placement(transformation(extent={{194,-62},{214,-42}}),
            iconTransformation(extent={{100,32},{120,52}})));
      Bodylight.Types.RealIO.PressureOutput pCO2 "CO2 partial pressure "
        annotation (Placement(transformation(extent={{194,-62},{214,-42}}),
            iconTransformation(extent={{100,-28},{120,-8}})));
      Bodylight.Types.RealIO.PressureOutput pN2 "N2 partial pressure "
        annotation (Placement(transformation(extent={{194,-62},{214,-42}}),
            iconTransformation(extent={{100,-86},{120,-66}})));
      Bodylight.Types.RealIO.PressureOutput pAr "Argon partial pressure "
        annotation (Placement(transformation(extent={{194,-62},{214,-42}}),
            iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={4,-110})));
      Bodylight.Types.RealIO.PressureOutput pH2O "wapor partial pressure "
        annotation (Placement(transformation(extent={{194,-62},{214,-42}}),
            iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-62,-110})));

      Modelica.Units.SI.Pressure pH2Os "saturated vapor pressure of water";
      Modelica.Units.SI.MolarVolume VmO2; //O2 molar volume
      Modelica.Units.SI.MolarVolume VmCO2; //CO2 molar volume
      Modelica.Units.SI.MolarVolume VmN2; //N2 molar volume
      Modelica.Units.SI.MolarVolume VmAr; //Argon molar volume
      Real tempC = temp-273.15; //temperature in Celsius
    equation
      //https://en.wikipedia.org/wiki/Arden_Buck_equation
      //reusult is in kPa
      //pH2O = 0.61121*exp((18.678 - tempC/234.5)*(tempC/(257.14 + tempC))) in kPa
      pH2Os = 0.61121*exp((18.678 - tempC/234.5)*(tempC/(257.14 + tempC))) *1000;
      pH2O=pH2Os*hum;
      FH2O=pH2O/PB;
      pO2=(PB-pH2O)*FO2d;
      pN2=(PB-pH2O)*FN2d;
      pCO2=(PB-pH2O)*FCO2d;
      pAr=(PB-pH2O)*FArd;
      FO2=pO2/PB;
      FCO2=pCO2/PB;
      FN2=pN2/PB;
      FAr=pAr/PB;

      VmO2stand*Pstand/273.15=VmO2*PB/temp;
      VmCO2stand*Pstand/273.15=VmCO2*PB/temp;
      VmN2stand*Pstand/273.15=VmN2*PB/temp;
      VmArstand*Pstand/273.15=VmAr*PB/temp;
      O2conc=FO2/VmO2;
      CO2conc=FCO2/VmCO2;
      N2conc=FN2/VmN2;
      Arconc=FAr/VmAr;
      annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
            Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={28,108,200},
              fillColor={85,255,255},
              fillPattern=FillPattern.Solid),
            Text(
              extent={{56,-50},{98,-84}},
              textColor={28,108,200},
              textString="N2"),
            Text(
              extent={{38,22},{92,-24}},
              textColor={28,108,200},
              textString="CO2"),
            Text(
              extent={{48,90},{94,54}},
              textColor={28,108,200},
              textString="O2"),
            Text(
              extent={{2,-64},{44,-96}},
              textColor={28,108,200},
              textString="Ar"),
            Text(
              extent={{-80,-58},{-26,-104}},
              textColor={28,108,200},
              textString="H2O"),
            Text(
              extent={{-80,-4},{-26,-50}},
              textColor={28,108,200},
              textString="hum"),
            Text(
              extent={{-76,56},{-22,10}},
              textColor={28,108,200},
              textString="temp"),
            Text(
              extent={{-78,102},{-22,54}},
              textColor={28,108,200},
              textString="press")}), Diagram(coordinateSystem(preserveAspectRatio=false)));
    end Athmosphere;

    model AlveolarVentilation
      Bodylight.Types.RealIO.VolumeFlowRateInput VAi annotation (Placement(
            transformation(extent={{106,22},{146,62}}), iconTransformation(extent={{
                -110,82},{-98,94}})));
      Bodylight.Types.RealIO.PressureInput pressure annotation (Placement(
            transformation(extent={{106,58},{146,98}}), iconTransformation(extent={{
                -110,62},{-96,76}})));
      Bodylight.Types.RealIO.TemperatureInput temperature annotation (Placement(
            transformation(extent={{104,2},{144,42}}), iconTransformation(extent={{-110,
                38},{-96,52}})));
      Bodylight.Types.RealIO.ConcentrationInput O2i_conc annotation (Placement(
            transformation(extent={{102,28},{142,68}}), iconTransformation(extent={{
                -110,14},{-96,28}})));
      Bodylight.Types.RealIO.ConcentrationInput CO2i_conc annotation (Placement(
            transformation(extent={{102,28},{142,68}}), iconTransformation(extent={{
                -112,-10},{-98,4}})));
      Bodylight.Types.Fraction FACO2;
      Bodylight.Types.Fraction FAO2;

      Bodylight.Types.MolarFlowRate O2i_flow;
      Bodylight.Types.MolarFlowRate CO2i_flow;
      Bodylight.Types.MolarFlowRate O2e_flow;
      Bodylight.Types.MolarFlowRate CO2e_flow;

      Bodylight.Types.VolumeFlowRate O2iSTPD_flow;
      Bodylight.Types.VolumeFlowRate CO2iSTPD_flow;
      Bodylight.Types.VolumeFlowRate O2eSTPD_flow;
      Bodylight.Types.VolumeFlowRate CO2eSTPD_flow;

      Bodylight.Types.VolumeFlowRate O2iBTPS_flow;
      Bodylight.Types.VolumeFlowRate CO2iBTPS_flow;
      Bodylight.Types.VolumeFlowRate O2eBTPS_flow;
      Bodylight.Types.VolumeFlowRate CO2eBTPS_flow;

      Bodylight.Types.RealIO.PressureOutput pAO2 annotation (Placement(
            transformation(extent={{132,40},{152,60}}), iconTransformation(extent={{
                102,68},{122,88}})));
      Bodylight.Types.RealIO.PressureOutput pACO2 annotation (Placement(
            transformation(extent={{132,40},{152,60}}), iconTransformation(extent={{102,46},
                {122,66}})));
      Bodylight.Types.RealIO.VolumeFlowRateOutput VAe annotation (Placement(
            transformation(extent={{148,-14},{168,6}}), iconTransformation(extent={{102,22},
                {122,42}})));

      Bodylight.Types.RealIO.PressureInput pH2O annotation (Placement(
            transformation(extent={{106,58},{146,98}}), iconTransformation(extent={{
                -112,-32},{-98,-18}})));

      Bodylight.Types.VolumeFlowRate VO2STPD;
      Bodylight.Types.VolumeFlowRate VCO2STPD;
      Bodylight.Types.VolumeFlowRate VAiSTPD;
      Bodylight.Types.VolumeFlowRate VO2iSTPD;
      Bodylight.Types.VolumeFlowRate VCO2iSTPD;
      Bodylight.Types.VolumeFlowRate VAeSTPD;
      Bodylight.Types.VolumeFlowRate VO2eSTPD;
      Bodylight.Types.VolumeFlowRate VCO2eSTPD;
      Bodylight.Types.Fraction FACO2dry;
      Bodylight.Types.Fraction FAO2dry;
      Bodylight.Types.Pressure PAO2;
      Bodylight.Types.Pressure PACO2;
      Bodylight.Types.VolumeFlowRate VAeBTPS;



      Bodylight.Types.VolumeFlowRate O2MSTPD;
      Bodylight.Types.VolumeFlowRate CO2MSTPD;

      parameter Modelica.Units.SI.MolarVolume VmO2stand=22.3922/1000
        "O2 molar volume at 0°C 760 mmHg";
      parameter Modelica.Units.SI.MolarVolume VmCO2stand=22.2629/1000
        "CO2 molar volume at 0°C 760 mmHg";
      parameter Modelica.Units.SI.Pressure Pstand=101325;
      //pressure 760 mmHg in Pa
      Real STPD_to_BTPS;

      parameter Bodylight.Types.Fraction FO2i_dry=0.2095;
      parameter Bodylight.Types.Fraction FCO2i_dry=0.04;

      Bodylight.Types.RealIO.MolarFlowRateInput O2M annotation (Placement(
            transformation(extent={{24,-52},{64,-12}}), iconTransformation(
            extent={{-12,-12},{12,12}},
            rotation=90,
            origin={-74,-106})));
      Bodylight.Types.RealIO.MolarFlowRateInput CO2M annotation (Placement(
            transformation(extent={{24,-52},{64,-12}}), iconTransformation(
            extent={{-12,-12},{12,12}},
            rotation=90,
            origin={46,-106})));
    equation
      //conversion O2M and CO2M from mmol/min to ml/min
      O2MSTPD=O2M*VmO2stand;
      CO2MSTPD=CO2M*VmCO2stand;


      //conversion coefficient from STPD to BTPS
      STPD_to_BTPS = Pstand*temperature/(273.15*(pressure - pH2O));

      //balance O2 and CO2 molar flows
      O2i_flow=VAi*O2i_conc;
      CO2i_flow=VAi*CO2i_conc;
      O2e_flow=O2i_flow-O2M;
      CO2e_flow=CO2i_flow+CO2M;

      //balance O2 and CO2 volume flows
      O2iSTPD_flow=O2i_flow*VmO2stand;
      CO2iSTPD_flow=CO2i_flow*VmCO2stand;
      O2eSTPD_flow=O2e_flow*VmO2stand;
      CO2eSTPD_flow=CO2e_flow*VmCO2stand;

      O2iBTPS_flow=O2iSTPD_flow*STPD_to_BTPS;
      CO2iBTPS_flow=CO2iSTPD_flow*STPD_to_BTPS;
      O2eBTPS_flow=O2eSTPD_flow*STPD_to_BTPS;
      CO2eBTPS_flow=CO2eSTPD_flow*STPD_to_BTPS;

      //Expired alveolar ventilation calculation
      //VAe=VAi+(O2eBTPS_flow-O2iBTPS_flow)+(CO2eBTPS_flow+CO2iBTPS_flow);
      VAe=VAi+(O2eBTPS_flow-O2iBTPS_flow)+(CO2eBTPS_flow-CO2iBTPS_flow);

      //Alveolar fraction in dry expired gas concentrations calculation
      FAO2=O2eSTPD_flow/(VAe/STPD_to_BTPS);
      FACO2=CO2eSTPD_flow/(VAe/STPD_to_BTPS);

      //Alveolat partial pressures calculation
      pACO2=FACO2*(pressure - pH2O);
      pAO2=FAO2*(pressure - pH2O);


      //***********************
      //alternative calculations
      //***********************

      //calculation of VO2 and VCO2 in STPD
      VO2STPD=O2M*VmO2stand;
      VCO2STPD=CO2M*VmCO2stand;

      //calculation of VAi_STPD
      VAiSTPD=VAi/STPD_to_BTPS;

      //Calculation of O2 and CO2 flows in inspired gas in STPD
      VO2iSTPD=VAiSTPD*FO2i_dry;
      VCO2iSTPD=VAiSTPD*FCO2i_dry;

      //Calculation of expired alveoar ventilation in STPD
      VAeSTPD=VAiSTPD-VO2STPD+VCO2STPD;

      //Calculation of expired O2 and CO2 flows in STPD
      VO2eSTPD=VO2iSTPD-VO2STPD;
      VCO2eSTPD=VCO2iSTPD+VCO2STPD;

      //Calculation of dry expired fractions of CO2 and O2
      FACO2dry=VCO2eSTPD/VAeSTPD;
      FAO2dry=VO2eSTPD/VAeSTPD;

      //Calsulation of alveolar PO2 and PCO2 (in BTPS conditions)
      PAO2=FAO2dry*(pressure-pH2O);
      PACO2=FACO2dry*(pressure-pH2O);

      //Calculation of expired BTPS aůveolar volume
      VAeBTPS=VAeSTPD*STPD_to_BTPS;

         annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
              Rectangle(
              extent={{-100,100},{102,-100}},
              lineColor={255,170,170},
              fillColor={0,255,255},
              fillPattern=FillPattern.Solid), Text(
              extent={{-72,66},{56,-30}},
              textColor={28,108,200},
              textString="AlvVentilation")}),   Diagram(coordinateSystem(
              preserveAspectRatio=false)));
    end AlveolarVentilation;

    model AlvVent
      Bodylight.Types.RealIO.VolumeFlowRateInput VAi_BTPS annotation (Placement(
            transformation(extent={{106,22},{146,62}}), iconTransformation(extent={{-110,76},
                {-98,88}})));
      Bodylight.Types.RealIO.PressureInput pressure annotation (Placement(
            transformation(extent={{106,58},{146,98}}), iconTransformation(extent={{-114,52},
                {-100,66}})));
      Bodylight.Types.RealIO.TemperatureInput temperature annotation (Placement(
            transformation(extent={{104,2},{144,42}}), iconTransformation(extent={{-114,26},
                {-100,40}})));


      Bodylight.Types.RealIO.PressureOutput pAO2 annotation (Placement(
            transformation(extent={{132,40},{152,60}}), iconTransformation(extent={{
                102,68},{122,88}})));
      Bodylight.Types.RealIO.PressureOutput pACO2 annotation (Placement(
            transformation(extent={{132,40},{152,60}}), iconTransformation(extent={{102,16},
                {122,36}})));
      Bodylight.Types.RealIO.VolumeFlowRateOutput VAe_BTPS annotation (Placement(
            transformation(extent={{148,-14},{168,6}}), iconTransformation(extent={{102,-36},
                {122,-16}})));

      Bodylight.Types.Pressure pH2O;

      Bodylight.Types.VolumeFlowRate VO2_STPD;
      Bodylight.Types.VolumeFlowRate VCO2_STPD;
      Bodylight.Types.VolumeFlowRate VAi_STPD;
      Bodylight.Types.VolumeFlowRate VO2i_STPD;
      Bodylight.Types.VolumeFlowRate VCO2i_STPD;
      Bodylight.Types.VolumeFlowRate VAe_STPD;
      Bodylight.Types.VolumeFlowRate VO2e_STPD;
      Bodylight.Types.VolumeFlowRate VCO2e_STPD;
      Bodylight.Types.Fraction FACO2dry;
      Bodylight.Types.Fraction FAO2dry;



      parameter Modelica.Units.SI.MolarVolume VmO2stand=22.3922/1000
        "O2 molar volume at 0°C 760 mmHg";
      parameter Modelica.Units.SI.MolarVolume VmCO2stand=22.2629/1000
        "CO2 molar volume at 0°C 760 mmHg";
      parameter Modelica.Units.SI.Pressure Pstand=101325;
      //pressure 760 mmHg in Pa
      Real STPD_to_BTPS;
      Real tempC = temperature-273.15; //temperature in Celsius

      Bodylight.Types.RealIO.MolarFlowRateInput O2M annotation (Placement(
            transformation(extent={{24,-52},{64,-12}}), iconTransformation(
            extent={{-9,-9},{9,9}},
            rotation=0,
            origin={-105,-59})));
      Bodylight.Types.RealIO.MolarFlowRateInput CO2M annotation (Placement(
            transformation(extent={{24,-52},{64,-12}}), iconTransformation(
            extent={{-9,-9},{9,9}},
            rotation=0,
            origin={-107,-91})));
      Bodylight.Types.RealIO.FractionInput FiO2d annotation (Placement(
            transformation(extent={{38,36},{78,76}}), iconTransformation(extent={{-116,-2},
                {-100,14}})));
      Bodylight.Types.RealIO.FractionInput FiCO2d annotation (Placement(
            transformation(extent={{38,36},{78,76}}), iconTransformation(extent={{-116,
                -36},{-100,-20}})));
    equation
      //conversion coefficient from STPD to BTPS
      STPD_to_BTPS = Pstand*temperature/(273.15*(pressure - pH2O));

      //https://en.wikipedia.org/wiki/Arden_Buck_equation
      //reusult is in kPa
      //pH2O = 0.61121*exp((18.678 - tempC/234.5)*(tempC/(257.14 + tempC))) in kPa

      pH2O = 0.61121*exp((18.678 - tempC/234.5)*(tempC/(257.14 + tempC))) *1000;

      //alternative calculations

      //calculation of O2 consumtion and CO2 production rates in STPD/s from mol/s
      VO2_STPD=O2M*VmO2stand;
      VCO2_STPD=CO2M*VmCO2stand;

      //calculation of VAi_STPD
      VAi_STPD=VAi_BTPS/STPD_to_BTPS;

      // calculation of inspired gases flow in l STPD/min
      VO2i_STPD=VAi_STPD*FiO2d;
      VCO2i_STPD=VAi_STPD*FiCO2d;

      //calculation of expired alveolar ventilation in l STPD/min
      VAe_STPD=VAi_STPD-VO2_STPD+VCO2_STPD;

      //calculation of expired gases flow in STPD
      VO2e_STPD=VO2i_STPD-VO2_STPD;
      VCO2e_STPD=VCO2i_STPD+VCO2_STPD;

      //calculation of FeO2 and FeCO2 in dry expired gas
      FACO2dry=VCO2e_STPD/VAe_STPD;
      FAO2dry=VO2e_STPD/VAe_STPD;

      //calculation of expired alveolar ventilation
      //in l BTPS/time (VAe_BTPS)
      VAe_BTPS=VAe_STPD*STPD_to_BTPS;

      ////calculation of PAO2 and PACO2 (at BTPS)
      pAO2=FAO2dry*(pressure-pH2O);
      pACO2=FACO2dry*(pressure-pH2O);

         annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
              Rectangle(
              extent={{-100,100},{102,-100}},
              lineColor={255,170,170},
              fillColor={0,255,255},
              fillPattern=FillPattern.Solid),
            Text(
              extent={{-96,98},{-54,66}},
              textColor={28,108,200},
              textString="VAiBTPS"),
            Text(
              extent={{-98,76},{-60,44}},
              textColor={28,108,200},
              textString="pressure"),
            Text(
              extent={{-96,52},{-40,14}},
              textColor={28,108,200},
              textString="temperature"),
            Text(
              extent={{-94,20},{-58,-8}},
              textColor={28,108,200},
              textString="FiO2d"),
            Text(
              extent={{-92,-12},{-56,-40}},
              textColor={28,108,200},
              textString="FiCO2d"),
            Text(
              extent={{-86,-50},{-50,-78}},
              textColor={28,108,200},
              textString="MO2"),
            Text(
              extent={{-88,-72},{-52,-100}},
              textColor={28,108,200},
              textString="MCO2"),
            Text(
              extent={{62,40},{98,12}},
              textColor={28,108,200},
              textString="PAO2"),
            Text(
              extent={{62,-12},{98,-40}},
              textColor={28,108,200},
              textString="PACO2"),
            Text(
              extent={{46,96},{90,64}},
              textColor={28,108,200},
              textString="VAeBTPS"),
            Text(
              extent={{-118,-128},{128,-172}},
              textColor={28,108,200},
              textString="%name")}),            Diagram(coordinateSystem(
              preserveAspectRatio=false)));
    end AlvVent;
  end Parts;

  package Tests
    model TestAthmosphere
      extends Modelica.Icons.Example;
      Parts.Athmosphere athmosphere
        annotation (Placement(transformation(extent={{-30,2},{6,38}})));
      Bodylight.Types.Constants.PressureConst pressure(k=101325.0144354)
        annotation (Placement(transformation(extent={{-78,28},{-70,36}})));
      Bodylight.Types.Constants.TemperatureConst temperature(k(displayUnit="degC")=
             310.15)
        annotation (Placement(transformation(extent={{-76,14},{-68,22}})));
      Bodylight.Types.Constants.FractionConst fraction(k=1)
        annotation (Placement(transformation(extent={{-72,0},{-64,8}})));
    equation
      connect(pressure.y, athmosphere.PB) annotation (Line(points={{-69,32},{-69,33.14},
              {-31.98,33.14}}, color={0,0,127}));
      connect(temperature.y, athmosphere.temp) annotation (Line(points={{-67,18},{-40,
              18},{-40,23.78},{-31.98,23.78}}, color={0,0,127}));
      connect(fraction.y, athmosphere.hum) annotation (Line(points={{-63,4},{-40,4},
              {-40,13.88},{-32.16,13.88}}, color={0,0,127}));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)));
    end TestAthmosphere;

    model TestALvVentilation
        extends Modelica.Icons.Example;
      Parts.Athmosphere athmosphere(Pstand(displayUnit="Torr"))
        annotation (Placement(transformation(extent={{-52,18},{-12,64}})));
      Bodylight.Types.Constants.PressureConst pressure(k(displayUnit="mmHg")
           = 101325.0144354)
        annotation (Placement(transformation(extent={{-96,54},{-88,62}})));
      Bodylight.Types.Constants.TemperatureConst temperature(k(displayUnit=
              "degC") = 310.15)
        annotation (Placement(transformation(extent={{-96,32},{-88,40}})));
      Bodylight.Types.Constants.FractionConst fraction(k(displayUnit="mmHg")
           = 1)
        annotation (Placement(transformation(extent={{-94,12},{-86,20}})));
      Bodylight.Types.Constants.VolumeFlowRateConst VAi(k=8.1e-05)
        annotation (Placement(transformation(extent={{-46,86},{-38,94}})));
      Bodylight.Types.Constants.MolarFlowRateConst MO2(k=
            0.00018666666666667)
        annotation (Placement(transformation(extent={{-12,-26},{2,-16}})));
      Bodylight.Types.Constants.MolarFlowRateConst MCO2(k=
            0.00016666666666667)
        annotation (Placement(transformation(extent={{-48,-60},{-34,-50}})));
      Parts.AlveolarVentilation alveolarVentilation(FCO2i_dry=0.0004)
        annotation (Placement(transformation(extent={{20,30},{76,84}})));
    equation
      connect(pressure.y, athmosphere.PB) annotation (Line(points={{-87,58},{
              -68,58},{-68,57.79},{-54.2,57.79}},
                               color={0,0,127}));
      connect(temperature.y, athmosphere.temp) annotation (Line(points={{-87,36},
              {-64,36},{-64,45.83},{-54.2,45.83}},
                                               color={0,0,127}));
      connect(fraction.y, athmosphere.hum) annotation (Line(points={{-85,16},{
              -62,16},{-62,33.18},{-54.4,33.18}},
                                           color={0,0,127}));

      connect(alveolarVentilation.VAi, VAi.y) annotation (Line(points={{18.88,
              80.76},{-37,80.76},{-37,90}},          color={0,0,127}));
      connect(alveolarVentilation.pressure, athmosphere.PB) annotation (Line(
            points={{19.16,75.63},{-64,75.63},{-64,57.79},{-54.2,57.79}}, color
            ={0,0,127}));
      connect(alveolarVentilation.temperature, athmosphere.temp) annotation (
          Line(points={{19.16,69.15},{-102,69.15},{-102,45.83},{-54.2,45.83}},
            color={0,0,127}));
      connect(athmosphere.O2conc, alveolarVentilation.O2i_conc) annotation (
          Line(points={{-10,55.26},{2,55.26},{2,62.67},{19.16,62.67}}, color={0,
              0,127}));
      connect(athmosphere.CO2conc, alveolarVentilation.CO2i_conc) annotation (
          Line(points={{-10,41.46},{4,41.46},{4,56.19},{18.6,56.19}}, color={0,
              0,127}));
      connect(athmosphere.pH2O, alveolarVentilation.pH2O) annotation (Line(
            points={{-44.4,15.7},{-44.4,4},{12,4},{12,50.25},{18.6,50.25}},
            color={0,0,127}));
      connect(MO2.y, alveolarVentilation.O2M) annotation (Line(points={{3.75,
              -21},{27.28,-21},{27.28,28.38}}, color={0,0,127}));
      connect(MCO2.y, alveolarVentilation.CO2M) annotation (Line(points={{
              -32.25,-55},{60.88,-55},{60.88,28.38}}, color={0,0,127}));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)));
    end TestALvVentilation;

    model TestAlvVent
          extends Modelica.Icons.Example
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)));
      Parts.AlvVent alvVent
        annotation (Placement(transformation(extent={{-10,-50},{50,0}})));
      Bodylight.Types.Constants.VolumeFlowRateConst VAi(k=8.1e-05)
        annotation (Placement(transformation(extent={{-78,26},{-70,34}})));
      Bodylight.Types.Constants.PressureConst pressure(k(displayUnit="mmHg") =
          101325.0144354)
        annotation (Placement(transformation(extent={{-68,-8},{-60,0}})));
      Bodylight.Types.Constants.TemperatureConst temperature(k(displayUnit=
              "degC") = 310.15)
        annotation (Placement(transformation(extent={{-46,-18},{-38,-10}})));
      Bodylight.Types.Constants.FractionConst FiO2d(k(displayUnit="%") = 0.2095)
        annotation (Placement(transformation(extent={{-78,-28},{-70,-20}})));
      Bodylight.Types.Constants.FractionConst FiCO2d(k(displayUnit="%") =
          0.0004)
        annotation (Placement(transformation(extent={{-56,-36},{-48,-28}})));
      Bodylight.Types.Constants.MolarFlowRateConst MCO2(k=0.00016666666666667)
        annotation (Placement(transformation(extent={{-74,-80},{-60,-70}})));
      Bodylight.Types.Constants.MolarFlowRateConst MO2(k=0.00018666666666667)
        annotation (Placement(transformation(extent={{-78,-58},{-64,-48}})));
      Modelica.Blocks.Sources.Ramp ramp(
        height=140,
        duration=140,
        offset=40)
        annotation (Placement(transformation(extent={{-88,56},{-68,76}})));
      Modelica.Blocks.Math.Product product1
        annotation (Placement(transformation(extent={{20,50},{40,70}})));
      Modelica.Blocks.Math.Gain gain(k(displayUnit="%") = 0.01)
        annotation (Placement(transformation(extent={{-48,56},{-28,76}})));
    equation
      connect(pressure.y, alvVent.pressure) annotation (Line(points={{-59,-4},{
              -50,-4},{-50,-6},{-20,-6},{-20,-10.25},{-12.1,-10.25}},
                                                  color={0,0,127}));
      connect(temperature.y, alvVent.temperature) annotation (Line(points={{-37,-14},
              {-20,-14},{-20,-16.75},{-12.1,-16.75}},  color={0,0,127}));
      connect(FiO2d.y, alvVent.FiO2d) annotation (Line(points={{-69,-24},{-40.7,
              -24},{-40.7,-23.5},{-12.4,-23.5}},
                                        color={0,0,127}));
      connect(FiCO2d.y, alvVent.FiCO2d)
        annotation (Line(points={{-47,-32},{-12.4,-32}},
                                                       color={0,0,127}));
      connect(MCO2.y, alvVent.CO2M) annotation (Line(points={{-58.25,-75},{
              -58.25,-76},{-20,-76},{-20,-47.75},{-12.1,-47.75}},
                                                                color={0,0,127}));
      connect(MO2.y, alvVent.O2M) annotation (Line(points={{-62.25,-53},{-62.25,
              -54},{-22,-54},{-22,-39.75},{-11.5,-39.75}},
                                          color={0,0,127}));
      connect(VAi.y, product1.u2) annotation (Line(points={{-69,30},{10,30},{10,
              54},{18,54}}, color={0,0,127}));
      connect(ramp.y, gain.u)
        annotation (Line(points={{-67,66},{-50,66}}, color={0,0,127}));
      connect(gain.y, product1.u1)
        annotation (Line(points={{-27,66},{18,66}}, color={0,0,127}));
      connect(alvVent.VAi_BTPS, product1.y) annotation (Line(points={{-11.2,
              -4.5},{-18,-4.5},{-18,8},{58,8},{58,60},{41,60}}, color={0,0,127}));
    end TestAlvVent;
  end Tests;
  annotation (uses(Modelica(version="4.0.0"), Bodylight(version="1.0")));
end RespirationModeling;
