classdef srv_app_matlab_raw < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                        matlab.ui.Figure
        InterfaceRecombinationCalculatorLabel  matlab.ui.control.Label
        SiButtonGroup                   matlab.ui.container.ButtonGroup
        nButton                         matlab.ui.control.RadioButton
        pButton                         matlab.ui.control.RadioButton
        UIAxes                          matlab.ui.control.UIAxes
        UIAxes2                         matlab.ui.control.UIAxes
        CalculateButton                 matlab.ui.control.Button
        XaxisDropDown_2                 matlab.ui.control.DropDown
        SweepButton                     matlab.ui.control.Button
        UsLabel                         matlab.ui.control.Label
        nLabel                          matlab.ui.control.Label
        pLabel                          matlab.ui.control.Label
        ExportButton                    matlab.ui.control.Button
        ClearButton                     matlab.ui.control.Button
        DielectricDropDownLabel         matlab.ui.control.Label
        DielectricDropDown              matlab.ui.control.DropDown
        Ndopant1cm3EditFieldLabel       matlab.ui.control.Label
        dopantField                     matlab.ui.control.NumericEditField
        n1cm3EditFieldLabel             matlab.ui.control.Label
        n1cm3EditField                  matlab.ui.control.NumericEditField
        Qfqcm2EditFieldLabel            matlab.ui.control.Label
        Qfqcm2EditField                 matlab.ui.control.NumericEditField
        XaxisDropDownLabel              matlab.ui.control.Label
        XaxisDropDown                   matlab.ui.control.DropDown
        minEditFieldLabel               matlab.ui.control.Label
        minEditField                    matlab.ui.control.NumericEditField
        maxEditFieldLabel               matlab.ui.control.Label
        maxEditField                    matlab.ui.control.NumericEditField
        YaxisDropDownLabel              matlab.ui.control.Label
        YaxisDropDown                   matlab.ui.control.DropDown
        MetalDropDownLabel              matlab.ui.control.Label
        MetalDropDown                   matlab.ui.control.DropDown
        qqcm2EditFieldLabel             matlab.ui.control.Label
        qqcm2EditField                  matlab.ui.control.NumericEditField
        maxEditField_2Label             matlab.ui.control.Label
        maxEditField_2                  matlab.ui.control.NumericEditField
        minEditField_2Label             matlab.ui.control.Label
        minEditField_2                  matlab.ui.control.NumericEditField
        Label_3                         matlab.ui.control.Label
        NumSweepEditField               matlab.ui.control.NumericEditField
        rEditField_2Label               matlab.ui.control.Label
        rEditField_2                    matlab.ui.control.NumericEditField
        dnmEditFieldLabel               matlab.ui.control.Label
        dnmEditField                    matlab.ui.control.NumericEditField
        E_fermieVEditFieldLabel         matlab.ui.control.Label
        E_fermieVEditField              matlab.ui.control.NumericEditField
        mEditField_2Label               matlab.ui.control.Label
        mEditField_2                    matlab.ui.control.NumericEditField
        mseVEditFieldLabel              matlab.ui.control.Label
        mseVEditField                   matlab.ui.control.NumericEditField
        seVEditFieldLabel               matlab.ui.control.Label
        seVEditField                    matlab.ui.control.NumericEditField
        WaferThicknesscmEditFieldLabel  matlab.ui.control.Label
        WaferThicknesscmEditField       matlab.ui.control.NumericEditField
        SeffcmsEditFieldLabel           matlab.ui.control.Label
        SeffcmsEditField                matlab.ui.control.NumericEditField
        Defects1to6Panel                matlab.ui.container.Panel
        typeDropDownLabel               matlab.ui.control.Label
        typeDropDown                    matlab.ui.control.DropDown
        DitEditFieldLabel               matlab.ui.control.Label
        DitEditField                    matlab.ui.control.NumericEditField
        EEtEditFieldLabel               matlab.ui.control.Label
        EEtEditField                    matlab.ui.control.NumericEditField
        DitLabel                        matlab.ui.control.Label
        DitEditField_2                  matlab.ui.control.NumericEditField
        nEditField_2Label               matlab.ui.control.Label
        nEditField_2                    matlab.ui.control.NumericEditField
        pEditFieldLabel                 matlab.ui.control.Label
        pEditField                      matlab.ui.control.NumericEditField
        typeDropDown_2Label             matlab.ui.control.Label
        typeDropDown_2                  matlab.ui.control.DropDown
        DitEditField_3Label             matlab.ui.control.Label
        DitEditField_3                  matlab.ui.control.NumericEditField
        EEtEditField_2Label             matlab.ui.control.Label
        EEtEditField_2                  matlab.ui.control.NumericEditField
        DitEditField_4Label             matlab.ui.control.Label
        DitEditField_4                  matlab.ui.control.NumericEditField
        nEditField_3Label               matlab.ui.control.Label
        nEditField_3                    matlab.ui.control.NumericEditField
        pEditField_2Label               matlab.ui.control.Label
        pEditField_2                    matlab.ui.control.NumericEditField
        typeDropDown_3Label             matlab.ui.control.Label
        typeDropDown_3                  matlab.ui.control.DropDown
        typeDropDown_4Label             matlab.ui.control.Label
        typeDropDown_4                  matlab.ui.control.DropDown
        typeDropDown_5Label             matlab.ui.control.Label
        typeDropDown_5                  matlab.ui.control.DropDown
        DitEditField_5Label             matlab.ui.control.Label
        DitEditField_5                  matlab.ui.control.NumericEditField
        EEtEditField_3Label             matlab.ui.control.Label
        EEtEditField_3                  matlab.ui.control.NumericEditField
        DitEditField_6Label             matlab.ui.control.Label
        DitEditField_6                  matlab.ui.control.NumericEditField
        nEditField_4Label               matlab.ui.control.Label
        nEditField_4                    matlab.ui.control.NumericEditField
        pEditField_3Label               matlab.ui.control.Label
        pEditField_3                    matlab.ui.control.NumericEditField
        DitEditField_7Label             matlab.ui.control.Label
        DitEditField_7                  matlab.ui.control.NumericEditField
        EEtEditField_4Label             matlab.ui.control.Label
        EEtEditField_4                  matlab.ui.control.NumericEditField
        DitEditField_8Label             matlab.ui.control.Label
        DitEditField_8                  matlab.ui.control.NumericEditField
        nEditField_5Label               matlab.ui.control.Label
        nEditField_5                    matlab.ui.control.NumericEditField
        pEditField_4Label               matlab.ui.control.Label
        pEditField_4                    matlab.ui.control.NumericEditField
        DitEditField_9Label             matlab.ui.control.Label
        DitEditField_9                  matlab.ui.control.NumericEditField
        EEtEditField_5Label             matlab.ui.control.Label
        EEtEditField_5                  matlab.ui.control.NumericEditField
        DitEditField_10Label            matlab.ui.control.Label
        DitEditField_10                 matlab.ui.control.NumericEditField
        nEditField_6Label               matlab.ui.control.Label
        nEditField_6                    matlab.ui.control.NumericEditField
        pEditField_5Label               matlab.ui.control.Label
        pEditField_5                    matlab.ui.control.NumericEditField
        typeDropDown_6Label             matlab.ui.control.Label
        typeDropDown_6                  matlab.ui.control.DropDown
        DitEditField_12Label            matlab.ui.control.Label
        DitEditField_12                 matlab.ui.control.NumericEditField
        nEditField_7Label               matlab.ui.control.Label
        nEditField_7                    matlab.ui.control.NumericEditField
        pEditField_6Label               matlab.ui.control.Label
        pEditField_6                    matlab.ui.control.NumericEditField
        EEtEditField_7Label             matlab.ui.control.Label
        EEtEditField_6                  matlab.ui.control.NumericEditField
        DitEditField_13Label            matlab.ui.control.Label
        DitEditField_11                 matlab.ui.control.NumericEditField
        ButtonGroup                     matlab.ui.container.ButtonGroup
        ChargedButton                   matlab.ui.control.RadioButton
        NeutralButton                   matlab.ui.control.RadioButton
        Label_4                         matlab.ui.control.Label
        Label_5                         matlab.ui.control.Label
        Label_6                         matlab.ui.control.Label
        Label_7                         matlab.ui.control.Label
        NumSweepEditField_2             matlab.ui.control.NumericEditField
        Sn0cmsEditFieldLabel            matlab.ui.control.Label
        Sn0cmsEditField                 matlab.ui.control.NumericEditField
        Sp0cmsEditFieldLabel            matlab.ui.control.Label
        Sp0cmsEditField                 matlab.ui.control.NumericEditField
        J0sfAcm2EditFieldLabel          matlab.ui.control.Label
        J0sfAcm2EditField               matlab.ui.control.NumericEditField
    end

    
    properties (Access = public)
        Ndop; %1/cm3 + for acceptor - for donors ---  1ohm cm=5e15
        tox; %oxide thicknes cm %%% 13.728e-6 for double layers
        DeltaWorkF; % metal to semicon work funciton diff [V]
        q=1.60e-19;%C electron charge
        epsrox=3.9;%insulator relative permitivity constant
        eps0=8.85e-14;%F/cm vacuum permitivity
        k=1.38e-23;%J/K boltazmann constant
        T=300;%K Temperature
        
        ni=1e10;%1/cm3 Semicon intrinsic carrier consentration
        phiM=4.1;%eV Metal work function
        epsrs=11.9;% Semicon relative permitivity
        chi=4.1; %V Semicon suceptibility 
        Eg=1.12; %eV Semicon gap energy
        m_e=9.1095e-31;%kg  eletron's mass
        h=6.626068e-34;% m2 kg / s plank's constant
        m_effC;% kg semicon electrons effective mass for conduction
        vth;
        phiF;
        Vt; % eV Thermal energy 
        n_0;
        p_0;
        Evector;
        Usrh;
        
        DitMatrix;
        CapXn;
        CapXp;
        DitAcc;
        DitDon;
        Vsurface;
        
        figAxes;
        
    end
    
    methods (Access = public)
        
        function results = recalculateall(app)
            
            if strcmp(app.SiButtonGroup.SelectedObject.Text,'n')
                app.Ndop=-app.dopantField.Value;
            else
                app.Ndop=app.dopantField.Value;
            end
            
            app.Vt=(app.k*app.T)/app.q;            
            app.m_effC=0.26*app.m_e;
            app.vth=sqrt(3*app.k*app.T/app.m_effC)/1e-2;%cm/s thermal velocity
            app.phiF=-(app.Ndop/abs(app.Ndop))*app.Vt*log(abs(app.Ndop/app.ni));%eV Semicon fermi level
            app.tox=app.dnmEditField.Value*1e-7;
            
            app.E_fermieVEditField.Value=app.phiF;
            
            
            WorkF=app.chi+app.Eg/2-app.phiF;
            app.seVEditField.Value=WorkF;
            
            switch app.MetalDropDown.Value
                case 'none'
                    app.DeltaWorkF=0-WorkF;
                otherwise 
                    app.DeltaWorkF=app.mEditField_2.Value-WorkF;
            end
            

            app.mseVEditField.Value=app.DeltaWorkF;
            
            if app.Ndop>0 %acceptors
                app.p_0=app.Ndop;%1/cm3 equilibrium holes  concentration ~N_a-
                app.n_0=app.ni^2/app.p_0;%1/cm3 equilibrium electrons concentration 
            else % donors
                app.n_0=-app.Ndop;%1/cm3 equilibrium holes  concentration ~N_d+
                app.p_0=app.ni^2/app.n_0;%1/cm3 equilibrium electrons concentration 
            end
            
           readoutDit(app)
     
           Cn=app.vth*sum(app.DitMatrix.*app.CapXn,1)';
           Cp=app.vth*sum(app.DitMatrix.*app.CapXp,1)';
           
           app.Sn0cmsEditField.Value=Cn(ceil(length(Cn)/2));
           app.Sp0cmsEditField.Value=Cp(ceil(length(Cp)/2));
           
    
           Delta_n=app.n1cm3EditField.Value;
           Qf=app.Qfqcm2EditField.Value;
            

           
           SigmaQ=app.qqcm2EditField.Value+1e7;

           [Recombination,app.Usrh,Phis,Charge,Q_si,Q_it,J0s]= Us0(app,Delta_n,Qf,SigmaQ,Cn,Cp);

           Seff=Recombination./Delta_n;
           app.SeffcmsEditField.Value=Seff;
           app.J0sfAcm2EditField.Value=J0s;

           replotall(app)
        end
        
        
        function readoutDit(app)
            
            app.Evector=linspace(-app.Eg/2,app.Eg/2,50)';
            app.DitAcc=zeros(size(app.Evector));
            app.DitDon=zeros(size(app.Evector));
            app.DitMatrix=zeros(6,length(app.Evector));
            app.CapXn=zeros(6,length(app.Evector))+1e-25;
            app.CapXp=zeros(6,length(app.Evector))+1e-25;
            
            DitParams=[app.DitEditField.Value;app.EEtEditField.Value;app.DitEditField_2.Value;app.nEditField_2.Value;app.pEditField.Value;app.Eg];
            app.DitMatrix(1,:)=calcDit(app,app.typeDropDown.Value,app.Evector,DitParams);
            app.CapXn(1,app.DitMatrix(1,:)>1e8)=DitParams(4);
            app.CapXp(1,app.DitMatrix(1,:)>1e8)=DitParams(5);
            
            DitParams=[app.DitEditField_3.Value;app.EEtEditField_2.Value;app.DitEditField_4.Value;app.nEditField_3.Value;app.pEditField_2.Value;app.Eg];
            app.DitMatrix(2,:)=calcDit(app,app.typeDropDown_2.Value,app.Evector,DitParams);
            app.CapXn(2,app.DitMatrix(2,:)>1e8)=DitParams(4);
            app.CapXp(2,app.DitMatrix(2,:)>1e8)=DitParams(5);

            
            DitParams=[app.DitEditField_5.Value;app.EEtEditField_3.Value;app.DitEditField_6.Value;app.nEditField_4.Value;app.pEditField_3.Value;app.Eg];
            app.DitMatrix(3,:)=calcDit(app,app.typeDropDown_3.Value,app.Evector,DitParams);
            app.CapXn(3,app.DitMatrix(3,:)>1e8)=DitParams(4);
            app.CapXp(3,app.DitMatrix(3,:)>1e8)=DitParams(5);
            
            DitParams=[app.DitEditField_7.Value;app.EEtEditField_4.Value;app.DitEditField_8.Value;app.nEditField_5.Value;app.pEditField_4.Value;app.Eg];
            app.DitMatrix(4,:)=calcDit(app,app.typeDropDown_4.Value,app.Evector,DitParams);
            app.CapXn(4,app.DitMatrix(4,:)>1e8)=DitParams(4);
            app.CapXp(4,app.DitMatrix(4,:)>1e8)=DitParams(5);
             
            DitParams=[app.DitEditField_9.Value;app.EEtEditField_5.Value;app.DitEditField_10.Value;app.nEditField_6.Value;app.pEditField_5.Value;app.Eg];
            app.DitMatrix(5,:)=calcDit(app,app.typeDropDown_5.Value,app.Evector,DitParams);
            app.CapXn(5,app.DitMatrix(5,:)>1e8)=DitParams(4);
            app.CapXp(5,app.DitMatrix(5,:)>1e8)=DitParams(5);
            
            DitParams=[app.DitEditField_11.Value;app.EEtEditField_6.Value;app.DitEditField_12.Value;app.nEditField_7.Value;app.pEditField_6.Value;app.Eg];
            app.DitMatrix(6,:)=calcDit(app,app.typeDropDown_6.Value,app.Evector,DitParams);
            app.CapXn(6,app.DitMatrix(6,:)>1e8)=DitParams(4);
            app.CapXp(6,app.DitMatrix(6,:)>1e8)=DitParams(5);
            
        end
        function DitVectorial = calcDit(app,distype, Evec,DitParams)

            DitTop=DitParams(1);
            Etrap=DitParams(2);
            DeltaDit=DitParams(3);
            capXn=DitParams(4);
            capXp=DitParams(5);
            Egap=DitParams(6);

        switch distype
            
            case 'Don -tail'
                
                m_v=(log(1e8)-log(DitTop))/(DeltaDit);
                DitVectorial=exp(((m_v)*Evec+(log(DitTop)-Etrap*(m_v))));
                app.DitDon=app.DitDon+DitVectorial;
                

            case 'Acc -tail'
               
                m_c=(log(1e8)-log(DitTop))/(-DeltaDit);
                DitVectorial=exp(((m_c)*Evec+(log(DitTop)-Etrap*(m_c))));
                app.DitAcc=app.DitAcc+DitVectorial;
                    
            case 'Acc -tophat'
                
                DitVectorial=DitTop*(heaviside(Evec-(Etrap-DeltaDit/2))-heaviside(Evec-(Etrap+DeltaDit/2)) )+1e-8;
                app.DitAcc=app.DitAcc+DitVectorial;
                
            case 'Don -tophat'
                
                DitVectorial=DitTop*(heaviside(Evec-(Etrap-DeltaDit/2))-heaviside(Evec-(Etrap+DeltaDit/2)) )+1e-8;
                app.DitDon=app.DitDon+DitVectorial;
                
            case 'Acc -gauss'
               
                DitVectorial=DitTop*exp(-(Evec-Etrap).^2/(2*DeltaDit^2));
                app.DitAcc=app.DitAcc+DitVectorial;
                
            case 'Don -gauss'
                
                DitVectorial=DitTop*exp(-(Evec-Etrap).^2/(2*DeltaDit^2));
                app.DitDon=app.DitDon+DitVectorial;
        
        end
        end
        
        function replotall(app)
            yyaxis(app.UIAxes,'left')

            
            app.UIAxes.set('DefaultLineLineWidth', 1.5);
            
        plot(app.UIAxes,app.Evector,app.DitMatrix(1,:));
        hold(app.UIAxes,'on')
        plot(app.UIAxes,app.Evector,app.DitMatrix(2,:));
        plot(app.UIAxes,app.Evector,app.DitMatrix(3,:));
        plot(app.UIAxes,app.Evector,app.DitMatrix(4,:));
        plot(app.UIAxes,app.Evector,app.DitMatrix(5,:));
        plot(app.UIAxes,app.Evector,app.DitMatrix(6,:));
        plot(app.UIAxes,app.Evector,app.Usrh,'b-');
        hold(app.UIAxes,'off')
        
         yyaxis(app.UIAxes,'right')
         app.UIAxes.set('DefaultLineMarkerSize', 8);
        plot(app.UIAxes,app.Evector,app.CapXn(1,:),'.');
        hold(app.UIAxes,'on')
        plot(app.UIAxes,app.Evector,app.CapXn(2,:),'.');
        plot(app.UIAxes,app.Evector,app.CapXn(3,:),'.');
        plot(app.UIAxes,app.Evector,app.CapXn(4,:),'.');
        plot(app.UIAxes,app.Evector,app.CapXn(5,:),'.');
        plot(app.UIAxes,app.Evector,app.CapXn(6,:),'.');
        
        plot(app.UIAxes,app.Evector,app.CapXp(1,:),'.m');
        plot(app.UIAxes,app.Evector,app.CapXp(2,:),'.m');
        plot(app.UIAxes,app.Evector,app.CapXp(3,:),'.m');
        plot(app.UIAxes,app.Evector,app.CapXp(4,:),'.m');
        plot(app.UIAxes,app.Evector,app.CapXp(5,:),'.m');
        plot(app.UIAxes,app.Evector,app.CapXp(6,:),'.m');
        hold(app.UIAxes,'off')
                    
        end
        
        function [F0,F1,F2,F3,F4,F5,F6]=Us0(app,Delta_n,Qf,SigmaQ,Cn,Cp)
            
            
            %excess concentration in stable region
            n_d=app.n_0+Delta_n;
            p_d=app.p_0+Delta_n;
            
            % cuasi fermi levels in stable region
            phi_n=-app.Vt*log(n_d/app.ni); 
            phi_p=app.Vt*log(p_d/app.ni); 
               
            %starting values
          
            E=app.Evector;
            
            n_1=app.ni*exp((E)/app.Vt); %1/cm3 SRH carrier concentrations at the interface
            p_1=app.ni*exp((-E)/app.Vt);
            
            Cp2=Cp;
            Cn2=Cn;
            Cr=Cp2./Cn2;
            
            if strcmp(app.ButtonGroup.SelectedObject.Text,'Charged')
                Dit_donors=app.DitDon;
                Dit_accpetors=app.DitAcc;
            else
                Dit_donors=0;Dit_accpetors=0;
            end
            
            Qsize=41;
            if Delta_n<2e13;Qsize=101;end

            Qf_vec=linspace(Qf-3*SigmaQ,Qf+3*SigmaQ,Qsize);
            PdistQ=(1./(SigmaQ*sqrt(2*pi))*exp(-(Qf_vec-Qf).^2./(2*SigmaQ^2)));
            
            PhisVar=zeros(Qsize,1);
            phiSran=[-(app.Eg/1.2+app.phiF),app.Eg/1.2-app.phiF];%initial range all bandgap

            for Qf_index=1:length(Qf_vec)
               
                Qerror=1e15;
                
                while min(Qerror)>1e4   
                    phi_s2=linspace(phiSran(1),phiSran(2));
                    Q_Si=-sign(phi_s2).*sqrt(2*app.k*app.T*app.ni*app.epsrs*app.eps0).*sqrt(exp((phi_p-phi_s2)./app.Vt)-exp((phi_p)/app.Vt)+ ...
                        exp((phi_s2-phi_n)./app.Vt)-exp((-phi_n)/app.Vt)+phi_s2.*(app.Ndop)/(app.Vt*app.ni))/app.q;
                    n_s=n_d*exp(phi_s2/app.Vt);p_s=p_d*exp(-phi_s2/app.Vt);
                    f_a=(n_s+Cr.*p_1)./((n_s+n_1)+Cr.*(p_s+p_1));
                    Qit_d=trapz(E,(1-f_a).*Dit_donors);
                    Qit_a=-trapz(E,(f_a).*Dit_accpetors);
                    Qerror=abs(Q_Si+Qit_a+Qit_d+Qf_vec(Qf_index));
                    [~,indPhi]=min(Qerror);
                    phiSran2=phiSran(2)-phiSran(1);
                    phiSran=[phi_s2(indPhi)-phiSran2/100,phi_s2(indPhi)+phiSran2/100];
                    if phiSran2<1e-8 && min(Qerror)>1e4
                        disp('possible convergence error')
                        break
                    end
                
                end
                
                phiSran=[phiSran(1)-SigmaQ*10e-13,phiSran(1)+SigmaQ*10e-13];% for next val of Qf
                
                phi_s=phi_s2(indPhi);
                PhisVar(Qf_index)=phi_s;

            end
            
            
            n_1=app.ni*exp(E./app.Vt); %1/cm3 SRH carrier concentrations at the interface
            p_1=app.ni*exp(-E./app.Vt);
            
            n_s=n_d*exp(PhisVar./app.Vt);
            p_s=p_d*exp(-PhisVar./app.Vt);
            U_SRH=zeros(length(PhisVar),1);
            for i=1:length(Qf_vec)
                U_SRH(i)= trapz(E,(n_s(i).*p_s(i)-app.ni^2).*((1)./((n_s(i)+n_1)./(Cp2)+(p_s(i)+p_1)./(Cn2))).*PdistQ(i));
                if i==25
                    U_srh2=(n_s(i).*p_s(i)-app.ni^2).*((1)./((n_s(i)+n_1)./(Cp2)+(p_s(i)+p_1)./(Cn2)));
                end
            end
            
            F0=trapz(Qf_vec,U_SRH);  
            F1=U_srh2;
            F2=PhisVar(ceil(Qsize/2));%surface potential in the semicon
            
            
            phi_s2=F2;
            Q_Si=-sign(phi_s2).*sqrt(2*app.k*app.T*app.ni*app.epsrs*app.eps0).*sqrt(exp((phi_p-phi_s2)./app.Vt)-exp((phi_p)/app.Vt)+ ...
                    exp((phi_s2-phi_n)./app.Vt)-exp((-phi_n)/app.Vt)+phi_s2.*(app.Ndop)/(app.Vt*app.ni))/app.q;
               
            n_s=n_d*exp(phi_s2/app.Vt);p_s=p_d*exp(-phi_s2/app.Vt);
            f_a=(n_s+Cr.*p_1)./((n_s+n_1)+Cr.*(p_s+p_1));
            Qit_d=trapz(E,(1-f_a).*Dit_donors);
            Qit_a=-trapz(E,(f_a).*Dit_accpetors);
            Qerror=abs(Q_Si+Qit_a+Qit_d+Qf_vec(Qf_index));
            
            F3=Qerror;%minimised charge
            F4=Q_Si;% charge on silicon
            F5=Qit_d+Qit_a;% charge on interface
            F6=app.q*F0./(((n_s.*p_s)./(app.ni.^2)-1))*1e15;%J0s in fA/cm2
        end
        
      
        
        
        
    end
    
    methods (Access = private)
        
        function [t_Aug,t_Rad]=IntrinsicTbulk(app,Ndop,Delta_n)
          
            %excess concentration in stable region
            n_d=app.n_0+Delta_n;
            p_d=app.p_0+Delta_n;
            
            
            Blow=4.73e-15 ; %[cm-3s-1] low injection radiative recombination probability
            % Values from Altermatt in NUSOD'2005
            Bmin=0.2+(0-0.2)/(1+(app.T/320).^(2.5));
            b1=1.5e18+(1e7-1.5e18)/(1+(app.T/550).^(3.0));
            b3=4e18+(1e9-4e18)/(1+(app.T/365).^(3.54));
            Brel=Bmin+(1.00-Bmin)./(1+((n_d+p_d)/b1).^(0.54)+((n_d+p_d)/b3).^(1.25));
            
            B=Brel.*Blow; %[cm-3s-1] radiative recombination probability
            Egnarrow=1e-3; % from J. Appl. Phys., Vol. 84, No. 7, 1 October 1998 Andreas Schenk
            nieff=app.ni*exp(Egnarrow/(2*app.k*app.T/app.q));
            geeh=1+13*(1-tanh((app.n_0/3.3e17).^(0.66)));
            gehh=1+7.5*(1-tanh((app.p_0/7e17).^(0.63)));
            t_intr=(Delta_n)./((n_d.*p_d-nieff.^2).*(2.5e-31*geeh*app.n_0+8.5e-32*gehh*app.p_0+3e-29*Delta_n.^(0.92)+B));
            t_Aug=(Delta_n)./((n_d.*p_d-nieff.^2).*(2.5e-31*geeh*app.n_0+8.5e-32*gehh*app.p_0+3e-29*Delta_n.^(0.92)));
            t_Rad=(Delta_n)./((n_d.*p_d-nieff.^2).*(B));
                      
        end
        
       
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
           app.UIAxes.FontSize=14;
           app.UIAxes2.FontSize=14;
           app.UIAxes.PlotBoxAspectRatioMode='auto';
           app.UIAxes2.PlotBoxAspectRatioMode='auto';
           
           
            app.UIAxes.YScale='log';
            app.UIAxes.YLim = [1e9 1e18];
            
            yyaxis(app.UIAxes,'right')
            app.UIAxes.YLabel.String='Capture Cross Sec (cm2)';
            app.UIAxes.YLim = [1e-18 1e-13];
            app.UIAxes.YScale='log';
            
            recalculateall(app)

        end

        % Callback function: DitEditField, DitEditField_10, 
        % DitEditField_11, DitEditField_12, DitEditField_2, 
        % DitEditField_3, DitEditField_4, DitEditField_5, 
        % DitEditField_6, DitEditField_7, DitEditField_8, 
        % DitEditField_9, EEtEditField, EEtEditField_2, 
        % EEtEditField_3, EEtEditField_4, EEtEditField_5, 
        % EEtEditField_6, SiButtonGroup, dopantField, mEditField_2, 
        % nEditField_2, nEditField_3, nEditField_4, nEditField_5, 
        % nEditField_6, nEditField_7, pEditField, pEditField_2, 
        % pEditField_3, pEditField_4, pEditField_5, pEditField_6, 
        % typeDropDown, typeDropDown_2, typeDropDown_3, 
        % typeDropDown_4, typeDropDown_5, typeDropDown_6
        function dopantFieldValueChanged(app, event)
            recalculateall(app)
        end

        % Button pushed function: CalculateButton, SweepButton
        function CalculateButtonPushed(app, event)
            doSweep=0;
            
            switch event.Source.Text                    
                case 'Sweep'
                    doSweep=1;
                    if strcmp(app.XaxisDropDown_2.Value,'Charge') && app.maxEditField_2.Value>1e13
                        app.maxEditField_2.Value=1e13;
                    end
                    SweepVar=logspace(log10(app.minEditField_2.Value),log10(app.maxEditField_2.Value),app.NumSweepEditField.Value)'; 
            end
            
            Cn=app.vth*sum(app.DitMatrix.*app.CapXn,1)';
            Cp=app.vth*sum(app.DitMatrix.*app.CapXp,1)';

            sizeSweep=app.NumSweepEditField_2.Value;
            
            switch app.XaxisDropDown.Value
                case 'Delta_n'   
                    Delta_n=logspace(log10(app.minEditField.Value),log10(app.maxEditField.Value),sizeSweep)';
                    Qf=app.Qfqcm2EditField.Value;
                    Xvar=Delta_n;
                    app.UIAxes2.XLabel.String='Excess carrier density (1/cm3)';app.UIAxes2.XScale='log';

                case 'Charge'
                    if app.maxEditField.Value>1e13
                        app.maxEditField.Value=1e13;
                    end
                    
                    Qf=[-logspace(log10(app.maxEditField.Value),log10(app.minEditField.Value),ceil(sizeSweep))';logspace(log10(app.minEditField.Value),log10(app.maxEditField.Value),ceil(sizeSweep))' ];
                    Delta_n=app.n1cm3EditField.Value;
                    Xvar=Qf-app.Qfqcm2EditField.Value;
                    app.UIAxes2.XLabel.String='External Charge (q/cm2)';app.UIAxes2.XScale='linear';
                    
                case 'Surface Voltage'
                    if app.maxEditField.Value>1e13
                        app.maxEditField.Value=1e13;
                    end
                    Qf=[-logspace(log10(app.maxEditField.Value),log10(app.minEditField.Value),ceil(sizeSweep))';logspace(log10(app.minEditField.Value),log10(app.maxEditField.Value),ceil(sizeSweep))' ];
                    Delta_n=app.n1cm3EditField.Value;
                    app.UIAxes2.XLabel.String='Dielectric surface voltage (V)';app.UIAxes2.XScale='linear';
            end
            
            if doSweep==1 && ~strcmp(app.XaxisDropDown.Value,app.XaxisDropDown_2.Value)
                switch app.XaxisDropDown_2.Value
                    case 'Delta_n' 
                        Delta_n=SweepVar;LegMatrix=cellstr(num2str(Delta_n,'%1.2e'));
                    case 'Charge'
                        Qf=[-flip(SweepVar);SweepVar];LegMatrix=cellstr(num2str(Qf,'%1.2e'));
                end
            else
                doSweep=0;
            end 
            
            
            Recombination=zeros(length(Delta_n),length(Qf));
            Phis=zeros(length(Delta_n),length(Qf));
            Charge=zeros(length(Delta_n),length(Qf));
            Seff=zeros(length(Delta_n),length(Qf));
            J0s=zeros(length(Delta_n),length(Qf));
            Q_it=zeros(length(Delta_n),length(Qf));
            Q_si=zeros(length(Delta_n),length(Qf));


            
            SigmaQ=app.qqcm2EditField.Value+1e8;
            % calculate the recombination rate
            for j=1:length(Qf)
                i=1;      
                for i=1:length(Delta_n)
                    [Recombination(i,j),~,Phis(i,j),Charge(i,j),Q_si(i,j),Q_it(i,j),J0s(i,j)]= Us0(app,Delta_n(i,1),Qf(j,1),SigmaQ,Cn,Cp);
                end
                Seff(:,j)=Recombination(:,j)./Delta_n;  
            end
            
            switch app.XaxisDropDown.Value
                case 'Surface Voltage'
                    Qoffset=app.Qfqcm2EditField.Value;
                    app.Vsurface=app.DeltaWorkF+Phis+(app.tox*(Qf-Qoffset)*app.q*ones(1,length(Delta_n))/(app.eps0*app.epsrox))';
                    Xvar=app.Vsurface;
            end
            
                    
            switch app.YaxisDropDown.Value
                case 'SRV'   
                    
                    plot(app.UIAxes2,Xvar,Seff);
                    hold(app.UIAxes2,'on')
                    
                    app.UIAxes2.YLabel.String='Seff (cm/s)';
                    app.UIAxes2.YScale='log';

                case 'Teff'
                    waferT=app.WaferThicknesscmEditField.Value; %[cm] wafer thickness
                    
                    [TauAuger,TauRad]=IntrinsicTbulk(app,app.Ndop,Delta_n);
                    
                    TauRadM=ones(size(Seff));
                    TauAugerM=ones(size(Seff));
                    for i=1:length(TauRad)
                       TauRadM(i,:)=TauRad(i);
                       TauAugerM(i,:)=TauAuger(i);
                    end
                    
                    TauEff=1./(1./(TauRadM)+1./TauAugerM+ 2*Seff./waferT);
                    
                    plot(app.UIAxes2,Xvar,TauEff);
                    
                    app.UIAxes2.YLabel.String='TauEff (s)';
                    app.UIAxes2.YScale='log';
                    
                case 'J0s'
                    plot(app.UIAxes2,Xvar,J0s);
                    hold(app.UIAxes2,'on')
                    
                    app.UIAxes2.YLabel.String='J0s (fA/cm^2)';
                    app.UIAxes2.YScale='log';

                    
            end
            hold(app.UIAxes2,'on')
            
%             app.UIAxes2.YMinorTick='on';
               if doSweep==1
                   
                   legend(app.UIAxes2,LegMatrix);
               end
            
        end

        % Value changed function: MetalDropDown
        function MetalDropDownValueChanged(app, event)
            value = app.MetalDropDown.Value;
            switch value
                case 'Al'
                    app.mEditField_2.Value=4.2;
                case 'Au'
                    app.mEditField_2.Value=5.1;
                case 'PEDOT'
                    app.mEditField_2.Value=5;
                case 'Ni'
                    app.mEditField_2.Value=5.2;
                case 'Ag'
                    app.mEditField_2.Value=4.15;
            end
            recalculateall(app);
        end

        % Value changed function: XaxisDropDown_2
        function XaxisDropDown_2ValueChanged(app, event)
            value = app.XaxisDropDown_2.Value;
            
            switch value
                case 'Charge'
                    app.maxEditField_2.Value=1e13;
                    app.minEditField_2.Value=1e9;
                case 'Delta_n'
                    app.maxEditField_2.Value=1e17;
                    app.minEditField_2.Value=1e13;
                    
            end
            
        end

        % Button pushed function: ExportButton
        function ExportButtonPushed(app, event)
            fignew = figure(1);
            
            if isempty(app.figAxes);app.figAxes = axes(fignew);end
            
            % Copy all UIAxes children, take over axes limits and aspect ratio.            
            allChildren = app.UIAxes2.XAxis.Parent.Children;
            copyobj(allChildren, app.figAxes)
            app.figAxes.XLim = app.UIAxes2.XLim;
            app.figAxes.YLim = app.UIAxes2.YLim;
            app.figAxes.XScale = app.UIAxes2.XScale;
            app.figAxes.YScale= app.UIAxes2.YScale;

            app.figAxes.DataAspectRatio = app.UIAxes2.DataAspectRatio;
            app.figAxes.NextPlot='add';
            
            xData=flip(get(app.UIAxes2.XAxis.Parent.Children,'Xdata')');
            yData=flip(get(app.UIAxes2.XAxis.Parent.Children,'Ydata')');
            if ~iscell(xData)
                assignin('base','xData',(xData));     
                assignin('base','yData',(yData));
            else
                assignin('base','xData',cell2mat(xData')');     
                assignin('base','yData',cell2mat(yData')');
            end
            open xData
            open yData

        end

        % Value changed function: XaxisDropDown, YaxisDropDown
        function YaxisDropDownValueChanged(app, event)
            value = app.YaxisDropDown.Value;
            
            switch app.XaxisDropDown.Value
                case 'Delta_n'
                    app.maxEditField.Value=1e17;
                    app.minEditField.Value=1e13;
                    
                otherwise
                    app.maxEditField.Value=1e13;
            end
            
        end

        % Button pushed function: ClearButton
        function ClearButtonPushed(app, event)
            hold(app.UIAxes2,'off');
            plot(app.UIAxes2,0);
        end

        % Value changed function: DielectricDropDown
        function DielectricDropDownValueChanged(app, event)
            value = app.DielectricDropDown.Value;
            
            switch value
                case 'SiO2'
                    app.rEditField_2.Value=3.9;
                    app.Qfqcm2EditField.Value=1e11;
                    
                case 'SiNx'
                    app.rEditField_2.Value=7.5;
                    app.Qfqcm2EditField.Value=1e12;
                    
                case 'AlOx'
                    app.rEditField_2.Value=8;
                    app.Qfqcm2EditField.Value=-1e12;
                    
                case 'a-Si:H'
                    app.rEditField_2.Value=11.5;
                    app.Qfqcm2EditField.Value=1e10;



            end
            recalculateall(app);
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 1053 691];
            app.UIFigure.Name = 'MATLAB App';

            % Create InterfaceRecombinationCalculatorLabel
            app.InterfaceRecombinationCalculatorLabel = uilabel(app.UIFigure);
            app.InterfaceRecombinationCalculatorLabel.FontSize = 26;
            app.InterfaceRecombinationCalculatorLabel.Position = [303 658 423 34];
            app.InterfaceRecombinationCalculatorLabel.Text = 'Interface Recombination Calculator';

            % Create SiButtonGroup
            app.SiButtonGroup = uibuttongroup(app.UIFigure);
            app.SiButtonGroup.SelectionChangedFcn = createCallbackFcn(app, @dopantFieldValueChanged, true);
            app.SiButtonGroup.Title = 'Si';
            app.SiButtonGroup.FontSize = 16;
            app.SiButtonGroup.Position = [11 484 100 45];

            % Create nButton
            app.nButton = uiradiobutton(app.SiButtonGroup);
            app.nButton.Text = 'n';
            app.nButton.Position = [8 1 29 22];
            app.nButton.Value = true;

            % Create pButton
            app.pButton = uiradiobutton(app.SiButtonGroup);
            app.pButton.Text = 'p';
            app.pButton.Position = [42 1 29 22];

            % Create UIAxes
            app.UIAxes = uiaxes(app.UIFigure);
            title(app.UIAxes, '')
            xlabel(app.UIAxes, 'Energy (eV)')
            ylabel(app.UIAxes, 'D_{it} (1/cm^{2}eV) and U_s(1/cm2.s.eV)')
            app.UIAxes.PlotBoxAspectRatio = [1 1.1702786377709 1];
            app.UIAxes.XLim = [-0.56 0.56];
            app.UIAxes.YLim = [100000000 1e+17];
            app.UIAxes.Box = 'on';
            app.UIAxes.YMinorTick = 'on';
            app.UIAxes.YScale = 'log';
            app.UIAxes.Position = [1 7 394 424];

            % Create UIAxes2
            app.UIAxes2 = uiaxes(app.UIFigure);
            title(app.UIAxes2, '')
            xlabel(app.UIAxes2, 'X')
            ylabel(app.UIAxes2, 'Y')
            app.UIAxes2.PlotBoxAspectRatio = [1.03767123287671 1 1];
            app.UIAxes2.Box = 'on';
            app.UIAxes2.YMinorTick = 'on';
            app.UIAxes2.Position = [542 7 411 363];

            % Create CalculateButton
            app.CalculateButton = uibutton(app.UIFigure, 'push');
            app.CalculateButton.ButtonPushedFcn = createCallbackFcn(app, @CalculateButtonPushed, true);
            app.CalculateButton.Position = [796 368 90 22];
            app.CalculateButton.Text = 'Calculate';

            % Create XaxisDropDown_2
            app.XaxisDropDown_2 = uidropdown(app.UIFigure);
            app.XaxisDropDown_2.Items = {'Delta_n', 'Charge'};
            app.XaxisDropDown_2.ValueChangedFcn = createCallbackFcn(app, @XaxisDropDown_2ValueChanged, true);
            app.XaxisDropDown_2.Position = [957 296 81 22];
            app.XaxisDropDown_2.Value = 'Delta_n';

            % Create SweepButton
            app.SweepButton = uibutton(app.UIFigure, 'push');
            app.SweepButton.ButtonPushedFcn = createCallbackFcn(app, @CalculateButtonPushed, true);
            app.SweepButton.Position = [994 177 48 22];
            app.SweepButton.Text = 'Sweep';

            % Create UsLabel
            app.UsLabel = uilabel(app.UIFigure);
            app.UsLabel.FontSize = 14;
            app.UsLabel.FontColor = [0 0 1];
            app.UsLabel.Position = [119 399 26 22];
            app.UsLabel.Text = 'Us ';

            % Create nLabel
            app.nLabel = uilabel(app.UIFigure);
            app.nLabel.FontSize = 16;
            app.nLabel.FontColor = [0.851 0.3255 0.098];
            app.nLabel.Position = [251 399 25 22];
            app.nLabel.Text = 'ÿn';

            % Create pLabel
            app.pLabel = uilabel(app.UIFigure);
            app.pLabel.FontSize = 16;
            app.pLabel.FontColor = [1 0 1];
            app.pLabel.Position = [301 399 25 22];
            app.pLabel.Text = 'ÿp';

            % Create ExportButton
            app.ExportButton = uibutton(app.UIFigure, 'push');
            app.ExportButton.ButtonPushedFcn = createCallbackFcn(app, @ExportButtonPushed, true);
            app.ExportButton.Position = [980 146 64 22];
            app.ExportButton.Text = 'Export';

            % Create ClearButton
            app.ClearButton = uibutton(app.UIFigure, 'push');
            app.ClearButton.ButtonPushedFcn = createCallbackFcn(app, @ClearButtonPushed, true);
            app.ClearButton.Position = [905 368 49 22];
            app.ClearButton.Text = 'Clear';

            % Create DielectricDropDownLabel
            app.DielectricDropDownLabel = uilabel(app.UIFigure);
            app.DielectricDropDownLabel.HorizontalAlignment = 'center';
            app.DielectricDropDownLabel.Position = [7 582 55 22];
            app.DielectricDropDownLabel.Text = 'Dielectric';

            % Create DielectricDropDown
            app.DielectricDropDown = uidropdown(app.UIFigure);
            app.DielectricDropDown.Items = {'SiO2', 'SiNx', 'AlOx', 'a-Si:H'};
            app.DielectricDropDown.ValueChangedFcn = createCallbackFcn(app, @DielectricDropDownValueChanged, true);
            app.DielectricDropDown.Position = [63 582 84 22];
            app.DielectricDropDown.Value = 'SiO2';

            % Create Ndopant1cm3EditFieldLabel
            app.Ndopant1cm3EditFieldLabel = uilabel(app.UIFigure);
            app.Ndopant1cm3EditFieldLabel.HorizontalAlignment = 'right';
            app.Ndopant1cm3EditFieldLabel.Position = [137 509 96 22];
            app.Ndopant1cm3EditFieldLabel.Text = 'Ndopant (1/cm3)';

            % Create dopantField
            app.dopantField = uieditfield(app.UIFigure, 'numeric');
            app.dopantField.Limits = [1000000000000 1e+20];
            app.dopantField.ValueDisplayFormat = '%2.2e';
            app.dopantField.ValueChangedFcn = createCallbackFcn(app, @dopantFieldValueChanged, true);
            app.dopantField.Position = [237 511 63 22];
            app.dopantField.Value = 5e+15;

            % Create n1cm3EditFieldLabel
            app.n1cm3EditFieldLabel = uilabel(app.UIFigure);
            app.n1cm3EditFieldLabel.HorizontalAlignment = 'right';
            app.n1cm3EditFieldLabel.Position = [144 480 74 22];
            app.n1cm3EditFieldLabel.Text = 'ÿn  (1/cm3)=';

            % Create n1cm3EditField
            app.n1cm3EditField = uieditfield(app.UIFigure, 'numeric');
            app.n1cm3EditField.Limits = [1000000000000 1e+20];
            app.n1cm3EditField.ValueDisplayFormat = '%2.2e';
            app.n1cm3EditField.Position = [225 480 75 22];
            app.n1cm3EditField.Value = 1e+15;

            % Create Qfqcm2EditFieldLabel
            app.Qfqcm2EditFieldLabel = uilabel(app.UIFigure);
            app.Qfqcm2EditFieldLabel.HorizontalAlignment = 'right';
            app.Qfqcm2EditFieldLabel.Position = [49 550 66 22];
            app.Qfqcm2EditFieldLabel.Text = 'Qf (q/cm2):';

            % Create Qfqcm2EditField
            app.Qfqcm2EditField = uieditfield(app.UIFigure, 'numeric');
            app.Qfqcm2EditField.Limits = [100000000 10000000000000];
            app.Qfqcm2EditField.ValueDisplayFormat = '%1.1e';
            app.Qfqcm2EditField.Position = [118 550 59 22];
            app.Qfqcm2EditField.Value = 100000000000;

            % Create XaxisDropDownLabel
            app.XaxisDropDownLabel = uilabel(app.UIFigure);
            app.XaxisDropDownLabel.HorizontalAlignment = 'right';
            app.XaxisDropDownLabel.Position = [592 400 39 22];
            app.XaxisDropDownLabel.Text = 'X-axis';

            % Create XaxisDropDown
            app.XaxisDropDown = uidropdown(app.UIFigure);
            app.XaxisDropDown.Items = {'Delta_n', 'Charge', 'Surface Voltage'};
            app.XaxisDropDown.ValueChangedFcn = createCallbackFcn(app, @YaxisDropDownValueChanged, true);
            app.XaxisDropDown.Position = [561 377 100 22];
            app.XaxisDropDown.Value = 'Charge';

            % Create minEditFieldLabel
            app.minEditFieldLabel = uilabel(app.UIFigure);
            app.minEditFieldLabel.HorizontalAlignment = 'right';
            app.minEditFieldLabel.Position = [787 398 24 22];
            app.minEditFieldLabel.Text = 'min';

            % Create minEditField
            app.minEditField = uieditfield(app.UIFigure, 'numeric');
            app.minEditField.Limits = [10000000000 1e+18];
            app.minEditField.ValueDisplayFormat = '%1.1e';
            app.minEditField.Position = [815 398 56 22];
            app.minEditField.Value = 10000000000;

            % Create maxEditFieldLabel
            app.maxEditFieldLabel = uilabel(app.UIFigure);
            app.maxEditFieldLabel.HorizontalAlignment = 'right';
            app.maxEditFieldLabel.Position = [882 398 28 22];
            app.maxEditFieldLabel.Text = 'max';

            % Create maxEditField
            app.maxEditField = uieditfield(app.UIFigure, 'numeric');
            app.maxEditField.Limits = [10000000000 1e+18];
            app.maxEditField.ValueDisplayFormat = '%1.1e';
            app.maxEditField.Position = [915 398 56 22];
            app.maxEditField.Value = 10000000000000;

            % Create YaxisDropDownLabel
            app.YaxisDropDownLabel = uilabel(app.UIFigure);
            app.YaxisDropDownLabel.HorizontalAlignment = 'right';
            app.YaxisDropDownLabel.Position = [706 400 38 22];
            app.YaxisDropDownLabel.Text = 'Y-axis';

            % Create YaxisDropDown
            app.YaxisDropDown = uidropdown(app.UIFigure);
            app.YaxisDropDown.Items = {'SRV', 'Teff', 'J0s'};
            app.YaxisDropDown.ValueChangedFcn = createCallbackFcn(app, @YaxisDropDownValueChanged, true);
            app.YaxisDropDown.Position = [675 378 100 22];
            app.YaxisDropDown.Value = 'SRV';

            % Create MetalDropDownLabel
            app.MetalDropDownLabel = uilabel(app.UIFigure);
            app.MetalDropDownLabel.HorizontalAlignment = 'center';
            app.MetalDropDownLabel.Position = [21 629 35 22];
            app.MetalDropDownLabel.Text = 'Metal';

            % Create MetalDropDown
            app.MetalDropDown = uidropdown(app.UIFigure);
            app.MetalDropDown.Items = {'none', 'Al', 'Au', 'Ni', 'Ag', 'PEDOT:PSS'};
            app.MetalDropDown.ValueChangedFcn = createCallbackFcn(app, @MetalDropDownValueChanged, true);
            app.MetalDropDown.Position = [63 629 84 22];
            app.MetalDropDown.Value = 'none';

            % Create qqcm2EditFieldLabel
            app.qqcm2EditFieldLabel = uilabel(app.UIFigure);
            app.qqcm2EditFieldLabel.HorizontalAlignment = 'right';
            app.qqcm2EditFieldLabel.Position = [183 550 69 22];
            app.qqcm2EditFieldLabel.Text = 'ÿ_q (q/cm2)';

            % Create qqcm2EditField
            app.qqcm2EditField = uieditfield(app.UIFigure, 'numeric');
            app.qqcm2EditField.Position = [267 550 58 22];

            % Create maxEditField_2Label
            app.maxEditField_2Label = uilabel(app.UIFigure);
            app.maxEditField_2Label.HorizontalAlignment = 'right';
            app.maxEditField_2Label.Position = [949 205 28 22];
            app.maxEditField_2Label.Text = 'max';

            % Create maxEditField_2
            app.maxEditField_2 = uieditfield(app.UIFigure, 'numeric');
            app.maxEditField_2.Limits = [10000000000 1e+18];
            app.maxEditField_2.ValueDisplayFormat = '%1.1e';
            app.maxEditField_2.Position = [982 205 56 22];
            app.maxEditField_2.Value = 1e+17;

            % Create minEditField_2Label
            app.minEditField_2Label = uilabel(app.UIFigure);
            app.minEditField_2Label.HorizontalAlignment = 'right';
            app.minEditField_2Label.Position = [957 237 24 22];
            app.minEditField_2Label.Text = 'min';

            % Create minEditField_2
            app.minEditField_2 = uieditfield(app.UIFigure, 'numeric');
            app.minEditField_2.Limits = [1000000000 1e+18];
            app.minEditField_2.ValueDisplayFormat = '%1.1e';
            app.minEditField_2.Position = [982 237 56.0784313725489 22];
            app.minEditField_2.Value = 10000000000000;

            % Create Label_3
            app.Label_3 = uilabel(app.UIFigure);
            app.Label_3.HorizontalAlignment = 'right';
            app.Label_3.Position = [980 267 25 22];
            app.Label_3.Text = '#';

            % Create NumSweepEditField
            app.NumSweepEditField = uieditfield(app.UIFigure, 'numeric');
            app.NumSweepEditField.Limits = [1 50];
            app.NumSweepEditField.ValueDisplayFormat = '%.0f';
            app.NumSweepEditField.Position = [1010 267 25 22];
            app.NumSweepEditField.Value = 3;

            % Create rEditField_2Label
            app.rEditField_2Label = uilabel(app.UIFigure);
            app.rEditField_2Label.HorizontalAlignment = 'right';
            app.rEditField_2Label.Position = [160 578 25 22];
            app.rEditField_2Label.Text = 'ÿr';

            % Create rEditField_2
            app.rEditField_2 = uieditfield(app.UIFigure, 'numeric');
            app.rEditField_2.Limits = [1 15];
            app.rEditField_2.Position = [196 578 30 22];
            app.rEditField_2.Value = 3.9;

            % Create dnmEditFieldLabel
            app.dnmEditFieldLabel = uilabel(app.UIFigure);
            app.dnmEditFieldLabel.HorizontalAlignment = 'right';
            app.dnmEditFieldLabel.Position = [234 578 39 22];
            app.dnmEditFieldLabel.Text = 'd (nm)';

            % Create dnmEditField
            app.dnmEditField = uieditfield(app.UIFigure, 'numeric');
            app.dnmEditField.Limits = [1 1000];
            app.dnmEditField.Position = [281 578 33 22];
            app.dnmEditField.Value = 100;

            % Create E_fermieVEditFieldLabel
            app.E_fermieVEditFieldLabel = uilabel(app.UIFigure);
            app.E_fermieVEditFieldLabel.HorizontalAlignment = 'right';
            app.E_fermieVEditFieldLabel.Position = [396 317 69 22];
            app.E_fermieVEditFieldLabel.Text = 'E_fermi (eV)';

            % Create E_fermieVEditField
            app.E_fermieVEditField = uieditfield(app.UIFigure, 'numeric');
            app.E_fermieVEditField.ValueDisplayFormat = '%2.2f';
            app.E_fermieVEditField.Editable = 'off';
            app.E_fermieVEditField.Position = [471 317 64 22];

            % Create mEditField_2Label
            app.mEditField_2Label = uilabel(app.UIFigure);
            app.mEditField_2Label.HorizontalAlignment = 'right';
            app.mEditField_2Label.Position = [166 621 25 22];
            app.mEditField_2Label.Text = 'ÿm';

            % Create mEditField_2
            app.mEditField_2 = uieditfield(app.UIFigure, 'numeric');
            app.mEditField_2.ValueChangedFcn = createCallbackFcn(app, @dopantFieldValueChanged, true);
            app.mEditField_2.Position = [206 621 69 22];
            app.mEditField_2.Value = 4.1;

            % Create mseVEditFieldLabel
            app.mseVEditFieldLabel = uilabel(app.UIFigure);
            app.mseVEditFieldLabel.HorizontalAlignment = 'right';
            app.mseVEditFieldLabel.Position = [403 258 62 22];
            app.mseVEditFieldLabel.Text = 'ÿÿms (eV)';

            % Create mseVEditField
            app.mseVEditField = uieditfield(app.UIFigure, 'numeric');
            app.mseVEditField.ValueDisplayFormat = '%2.2f';
            app.mseVEditField.Editable = 'off';
            app.mseVEditField.Position = [471 258 64 22];

            % Create seVEditFieldLabel
            app.seVEditFieldLabel = uilabel(app.UIFigure);
            app.seVEditFieldLabel.HorizontalAlignment = 'right';
            app.seVEditFieldLabel.Position = [421 288 44 22];
            app.seVEditFieldLabel.Text = 'ÿs (eV)';

            % Create seVEditField
            app.seVEditField = uieditfield(app.UIFigure, 'numeric');
            app.seVEditField.ValueDisplayFormat = '%2.2f';
            app.seVEditField.Editable = 'off';
            app.seVEditField.Position = [471 288 64 22];

            % Create WaferThicknesscmEditFieldLabel
            app.WaferThicknesscmEditFieldLabel = uilabel(app.UIFigure);
            app.WaferThicknesscmEditFieldLabel.HorizontalAlignment = 'right';
            app.WaferThicknesscmEditFieldLabel.Position = [107 448 124 22];
            app.WaferThicknesscmEditFieldLabel.Text = 'Wafer Thickness (cm):';

            % Create WaferThicknesscmEditField
            app.WaferThicknesscmEditField = uieditfield(app.UIFigure, 'numeric');
            app.WaferThicknesscmEditField.Limits = [0.001 1];
            app.WaferThicknesscmEditField.Position = [238 448 62 22];
            app.WaferThicknesscmEditField.Value = 0.02;

            % Create SeffcmsEditFieldLabel
            app.SeffcmsEditFieldLabel = uilabel(app.UIFigure);
            app.SeffcmsEditFieldLabel.HorizontalAlignment = 'right';
            app.SeffcmsEditFieldLabel.Position = [402 226 63 22];
            app.SeffcmsEditFieldLabel.Text = 'Seff (cm/s)';

            % Create SeffcmsEditField
            app.SeffcmsEditField = uieditfield(app.UIFigure, 'numeric');
            app.SeffcmsEditField.ValueDisplayFormat = '%2.2f';
            app.SeffcmsEditField.Editable = 'off';
            app.SeffcmsEditField.Position = [471 226 64 22];

            % Create Defects1to6Panel
            app.Defects1to6Panel = uipanel(app.UIFigure);
            app.Defects1to6Panel.TitlePosition = 'centertop';
            app.Defects1to6Panel.Title = 'Defects 1 to 6';
            app.Defects1to6Panel.Position = [340 435 682 224];

            % Create typeDropDownLabel
            app.typeDropDownLabel = uilabel(app.Defects1to6Panel);
            app.typeDropDownLabel.HorizontalAlignment = 'right';
            app.typeDropDownLabel.Visible = 'off';
            app.typeDropDownLabel.Position = [3 179 29 22];
            app.typeDropDownLabel.Text = 'type';

            % Create typeDropDown
            app.typeDropDown = uidropdown(app.Defects1to6Panel);
            app.typeDropDown.Items = {'Don -tail', 'Acc -tail', 'Acc -tophat', 'Don -tophat', 'Acc -gauss', 'Don -gauss'};
            app.typeDropDown.ValueChangedFcn = createCallbackFcn(app, @dopantFieldValueChanged, true);
            app.typeDropDown.Position = [3 179 99 22];
            app.typeDropDown.Value = 'Don -tail';

            % Create DitEditFieldLabel
            app.DitEditFieldLabel = uilabel(app.Defects1to6Panel);
            app.DitEditFieldLabel.HorizontalAlignment = 'right';
            app.DitEditFieldLabel.Position = [29 156 18 22];
            app.DitEditFieldLabel.Text = 'Dit';

            % Create DitEditField
            app.DitEditField = uieditfield(app.Defects1to6Panel, 'numeric');
            app.DitEditField.Limits = [100000000 1e+15];
            app.DitEditField.ValueChangedFcn = createCallbackFcn(app, @dopantFieldValueChanged, true);
            app.DitEditField.Position = [49 156 52.5333333333334 22];
            app.DitEditField.Value = 1e+15;

            % Create EEtEditFieldLabel
            app.EEtEditFieldLabel = uilabel(app.Defects1to6Panel);
            app.EEtEditFieldLabel.HorizontalAlignment = 'right';
            app.EEtEditFieldLabel.Position = [18 126 29 22];
            app.EEtEditFieldLabel.Text = 'E-Et';

            % Create EEtEditField
            app.EEtEditField = uieditfield(app.Defects1to6Panel, 'numeric');
            app.EEtEditField.Limits = [-1 1];
            app.EEtEditField.ValueChangedFcn = createCallbackFcn(app, @dopantFieldValueChanged, true);
            app.EEtEditField.Position = [49 126 52.5333333333334 22];
            app.EEtEditField.Value = -0.56;

            % Create DitLabel
            app.DitLabel = uilabel(app.Defects1to6Panel);
            app.DitLabel.HorizontalAlignment = 'right';
            app.DitLabel.Position = [19 94 28 22];
            app.DitLabel.Text = 'ÿDit';

            % Create DitEditField_2
            app.DitEditField_2 = uieditfield(app.Defects1to6Panel, 'numeric');
            app.DitEditField_2.Limits = [-2 2];
            app.DitEditField_2.ValueChangedFcn = createCallbackFcn(app, @dopantFieldValueChanged, true);
            app.DitEditField_2.Position = [49 94 52.5333333333334 22];
            app.DitEditField_2.Value = 0.12;

            % Create nEditField_2Label
            app.nEditField_2Label = uilabel(app.Defects1to6Panel);
            app.nEditField_2Label.HorizontalAlignment = 'right';
            app.nEditField_2Label.Position = [22 62 25 22];
            app.nEditField_2Label.Text = 'ÿ_n';

            % Create nEditField_2
            app.nEditField_2 = uieditfield(app.Defects1to6Panel, 'numeric');
            app.nEditField_2.Limits = [1e-20 1e-10];
            app.nEditField_2.ValueChangedFcn = createCallbackFcn(app, @dopantFieldValueChanged, true);
            app.nEditField_2.Position = [49 62 52.5333333333334 22];
            app.nEditField_2.Value = 2e-18;

            % Create pEditFieldLabel
            app.pEditFieldLabel = uilabel(app.Defects1to6Panel);
            app.pEditFieldLabel.HorizontalAlignment = 'right';
            app.pEditFieldLabel.Position = [22 30 25 22];
            app.pEditFieldLabel.Text = 'ÿ_p';

            % Create pEditField
            app.pEditField = uieditfield(app.Defects1to6Panel, 'numeric');
            app.pEditField.Limits = [1e-20 1e-10];
            app.pEditField.ValueChangedFcn = createCallbackFcn(app, @dopantFieldValueChanged, true);
            app.pEditField.Position = [49 30 52.5333333333334 22];
            app.pEditField.Value = 2e-18;

            % Create typeDropDown_2Label
            app.typeDropDown_2Label = uilabel(app.Defects1to6Panel);
            app.typeDropDown_2Label.HorizontalAlignment = 'right';
            app.typeDropDown_2Label.Visible = 'off';
            app.typeDropDown_2Label.Position = [109 176 29 22];
            app.typeDropDown_2Label.Text = 'type';

            % Create typeDropDown_2
            app.typeDropDown_2 = uidropdown(app.Defects1to6Panel);
            app.typeDropDown_2.Items = {'Don -tail', 'Acc -tail', 'Acc -tophat', 'Don -tophat', 'Acc -gauss', 'Don -gauss'};
            app.typeDropDown_2.ValueChangedFcn = createCallbackFcn(app, @dopantFieldValueChanged, true);
            app.typeDropDown_2.Position = [109 176 99 22];
            app.typeDropDown_2.Value = 'Don -tophat';

            % Create DitEditField_3Label
            app.DitEditField_3Label = uilabel(app.Defects1to6Panel);
            app.DitEditField_3Label.HorizontalAlignment = 'right';
            app.DitEditField_3Label.Position = [128 153 18 22];
            app.DitEditField_3Label.Text = 'Dit';

            % Create DitEditField_3
            app.DitEditField_3 = uieditfield(app.Defects1to6Panel, 'numeric');
            app.DitEditField_3.Limits = [100000000 1e+15];
            app.DitEditField_3.ValueChangedFcn = createCallbackFcn(app, @dopantFieldValueChanged, true);
            app.DitEditField_3.Position = [148 153 53 22];
            app.DitEditField_3.Value = 10000000000;

            % Create EEtEditField_2Label
            app.EEtEditField_2Label = uilabel(app.Defects1to6Panel);
            app.EEtEditField_2Label.HorizontalAlignment = 'right';
            app.EEtEditField_2Label.Position = [117 123 29 22];
            app.EEtEditField_2Label.Text = 'E-Et';

            % Create EEtEditField_2
            app.EEtEditField_2 = uieditfield(app.Defects1to6Panel, 'numeric');
            app.EEtEditField_2.Limits = [-1 1];
            app.EEtEditField_2.ValueChangedFcn = createCallbackFcn(app, @dopantFieldValueChanged, true);
            app.EEtEditField_2.Position = [148 123 53 22];
            app.EEtEditField_2.Value = -0.25;

            % Create DitEditField_4Label
            app.DitEditField_4Label = uilabel(app.Defects1to6Panel);
            app.DitEditField_4Label.HorizontalAlignment = 'right';
            app.DitEditField_4Label.Position = [118 91 28 22];
            app.DitEditField_4Label.Text = 'ÿDit';

            % Create DitEditField_4
            app.DitEditField_4 = uieditfield(app.Defects1to6Panel, 'numeric');
            app.DitEditField_4.Limits = [-2 2];
            app.DitEditField_4.ValueChangedFcn = createCallbackFcn(app, @dopantFieldValueChanged, true);
            app.DitEditField_4.Position = [148 91 53 22];
            app.DitEditField_4.Value = 0.5;

            % Create nEditField_3Label
            app.nEditField_3Label = uilabel(app.Defects1to6Panel);
            app.nEditField_3Label.HorizontalAlignment = 'right';
            app.nEditField_3Label.Position = [121 59 25 22];
            app.nEditField_3Label.Text = 'ÿ_n';

            % Create nEditField_3
            app.nEditField_3 = uieditfield(app.Defects1to6Panel, 'numeric');
            app.nEditField_3.Limits = [1e-20 1e-10];
            app.nEditField_3.ValueChangedFcn = createCallbackFcn(app, @dopantFieldValueChanged, true);
            app.nEditField_3.Position = [148 59 53 22];
            app.nEditField_3.Value = 1e-15;

            % Create pEditField_2Label
            app.pEditField_2Label = uilabel(app.Defects1to6Panel);
            app.pEditField_2Label.HorizontalAlignment = 'right';
            app.pEditField_2Label.Position = [121 27 25 22];
            app.pEditField_2Label.Text = 'ÿ_p';

            % Create pEditField_2
            app.pEditField_2 = uieditfield(app.Defects1to6Panel, 'numeric');
            app.pEditField_2.Limits = [1e-20 1e-10];
            app.pEditField_2.ValueChangedFcn = createCallbackFcn(app, @dopantFieldValueChanged, true);
            app.pEditField_2.Position = [148 27 53 22];
            app.pEditField_2.Value = 1e-16;

            % Create typeDropDown_3Label
            app.typeDropDown_3Label = uilabel(app.Defects1to6Panel);
            app.typeDropDown_3Label.HorizontalAlignment = 'right';
            app.typeDropDown_3Label.Visible = 'off';
            app.typeDropDown_3Label.Position = [217 176 29 22];
            app.typeDropDown_3Label.Text = 'type';

            % Create typeDropDown_3
            app.typeDropDown_3 = uidropdown(app.Defects1to6Panel);
            app.typeDropDown_3.Items = {'Don -tail', 'Acc -tail', 'Acc -tophat', 'Don -tophat', 'Acc -gauss', 'Don -gauss'};
            app.typeDropDown_3.ValueChangedFcn = createCallbackFcn(app, @dopantFieldValueChanged, true);
            app.typeDropDown_3.Position = [217 176 99 22];
            app.typeDropDown_3.Value = 'Acc -tophat';

            % Create typeDropDown_4Label
            app.typeDropDown_4Label = uilabel(app.Defects1to6Panel);
            app.typeDropDown_4Label.HorizontalAlignment = 'right';
            app.typeDropDown_4Label.Visible = 'off';
            app.typeDropDown_4Label.Position = [332 176 29 22];
            app.typeDropDown_4Label.Text = 'type';

            % Create typeDropDown_4
            app.typeDropDown_4 = uidropdown(app.Defects1to6Panel);
            app.typeDropDown_4.Items = {'Don -tail', 'Acc -tail', 'Acc -tophat', 'Don -tophat', 'Acc -gauss', 'Don -gauss'};
            app.typeDropDown_4.ValueChangedFcn = createCallbackFcn(app, @dopantFieldValueChanged, true);
            app.typeDropDown_4.Position = [332 176 99 22];
            app.typeDropDown_4.Value = 'Acc -tail';

            % Create typeDropDown_5Label
            app.typeDropDown_5Label = uilabel(app.Defects1to6Panel);
            app.typeDropDown_5Label.HorizontalAlignment = 'right';
            app.typeDropDown_5Label.Visible = 'off';
            app.typeDropDown_5Label.Position = [453 176 29 22];
            app.typeDropDown_5Label.Text = 'type';

            % Create typeDropDown_5
            app.typeDropDown_5 = uidropdown(app.Defects1to6Panel);
            app.typeDropDown_5.Items = {'Don -tail', 'Acc -tail', 'Acc -tophat', 'Don -tophat', 'Acc -gauss', 'Don -gauss'};
            app.typeDropDown_5.ValueChangedFcn = createCallbackFcn(app, @dopantFieldValueChanged, true);
            app.typeDropDown_5.Position = [453 176 99 22];
            app.typeDropDown_5.Value = 'Acc -tophat';

            % Create DitEditField_5Label
            app.DitEditField_5Label = uilabel(app.Defects1to6Panel);
            app.DitEditField_5Label.HorizontalAlignment = 'right';
            app.DitEditField_5Label.Position = [230 153 18 22];
            app.DitEditField_5Label.Text = 'Dit';

            % Create DitEditField_5
            app.DitEditField_5 = uieditfield(app.Defects1to6Panel, 'numeric');
            app.DitEditField_5.Limits = [100000000 1e+15];
            app.DitEditField_5.ValueChangedFcn = createCallbackFcn(app, @dopantFieldValueChanged, true);
            app.DitEditField_5.Position = [250 153 53 22];
            app.DitEditField_5.Value = 10000000000;

            % Create EEtEditField_3Label
            app.EEtEditField_3Label = uilabel(app.Defects1to6Panel);
            app.EEtEditField_3Label.HorizontalAlignment = 'right';
            app.EEtEditField_3Label.Position = [219 123 29 22];
            app.EEtEditField_3Label.Text = 'E-Et';

            % Create EEtEditField_3
            app.EEtEditField_3 = uieditfield(app.Defects1to6Panel, 'numeric');
            app.EEtEditField_3.Limits = [-1 1];
            app.EEtEditField_3.ValueChangedFcn = createCallbackFcn(app, @dopantFieldValueChanged, true);
            app.EEtEditField_3.Position = [250 123 53 22];
            app.EEtEditField_3.Value = 0.25;

            % Create DitEditField_6Label
            app.DitEditField_6Label = uilabel(app.Defects1to6Panel);
            app.DitEditField_6Label.HorizontalAlignment = 'right';
            app.DitEditField_6Label.Position = [220 91 28 22];
            app.DitEditField_6Label.Text = 'ÿDit';

            % Create DitEditField_6
            app.DitEditField_6 = uieditfield(app.Defects1to6Panel, 'numeric');
            app.DitEditField_6.Limits = [-2 2];
            app.DitEditField_6.ValueChangedFcn = createCallbackFcn(app, @dopantFieldValueChanged, true);
            app.DitEditField_6.Position = [250 91 53 22];
            app.DitEditField_6.Value = 0.5;

            % Create nEditField_4Label
            app.nEditField_4Label = uilabel(app.Defects1to6Panel);
            app.nEditField_4Label.HorizontalAlignment = 'right';
            app.nEditField_4Label.Position = [223 59 25 22];
            app.nEditField_4Label.Text = 'ÿ_n';

            % Create nEditField_4
            app.nEditField_4 = uieditfield(app.Defects1to6Panel, 'numeric');
            app.nEditField_4.Limits = [1e-20 1e-10];
            app.nEditField_4.ValueChangedFcn = createCallbackFcn(app, @dopantFieldValueChanged, true);
            app.nEditField_4.Position = [250 59 53 22];
            app.nEditField_4.Value = 1e-15;

            % Create pEditField_3Label
            app.pEditField_3Label = uilabel(app.Defects1to6Panel);
            app.pEditField_3Label.HorizontalAlignment = 'right';
            app.pEditField_3Label.Position = [223 27 25 22];
            app.pEditField_3Label.Text = 'ÿ_p';

            % Create pEditField_3
            app.pEditField_3 = uieditfield(app.Defects1to6Panel, 'numeric');
            app.pEditField_3.Limits = [1e-20 1e-10];
            app.pEditField_3.ValueChangedFcn = createCallbackFcn(app, @dopantFieldValueChanged, true);
            app.pEditField_3.Position = [250 27 53 22];
            app.pEditField_3.Value = 1e-16;

            % Create DitEditField_7Label
            app.DitEditField_7Label = uilabel(app.Defects1to6Panel);
            app.DitEditField_7Label.HorizontalAlignment = 'right';
            app.DitEditField_7Label.Position = [339 153 18 22];
            app.DitEditField_7Label.Text = 'Dit';

            % Create DitEditField_7
            app.DitEditField_7 = uieditfield(app.Defects1to6Panel, 'numeric');
            app.DitEditField_7.Limits = [100000000 1e+15];
            app.DitEditField_7.ValueChangedFcn = createCallbackFcn(app, @dopantFieldValueChanged, true);
            app.DitEditField_7.Position = [359 153 53 22];
            app.DitEditField_7.Value = 1e+15;

            % Create EEtEditField_4Label
            app.EEtEditField_4Label = uilabel(app.Defects1to6Panel);
            app.EEtEditField_4Label.HorizontalAlignment = 'right';
            app.EEtEditField_4Label.Position = [328 123 29 22];
            app.EEtEditField_4Label.Text = 'E-Et';

            % Create EEtEditField_4
            app.EEtEditField_4 = uieditfield(app.Defects1to6Panel, 'numeric');
            app.EEtEditField_4.Limits = [-1 1];
            app.EEtEditField_4.ValueChangedFcn = createCallbackFcn(app, @dopantFieldValueChanged, true);
            app.EEtEditField_4.Position = [359 123 53 22];
            app.EEtEditField_4.Value = 0.56;

            % Create DitEditField_8Label
            app.DitEditField_8Label = uilabel(app.Defects1to6Panel);
            app.DitEditField_8Label.HorizontalAlignment = 'right';
            app.DitEditField_8Label.Position = [329 91 28 22];
            app.DitEditField_8Label.Text = 'ÿDit';

            % Create DitEditField_8
            app.DitEditField_8 = uieditfield(app.Defects1to6Panel, 'numeric');
            app.DitEditField_8.Limits = [-2 2];
            app.DitEditField_8.ValueChangedFcn = createCallbackFcn(app, @dopantFieldValueChanged, true);
            app.DitEditField_8.Position = [359 91 53 22];
            app.DitEditField_8.Value = 0.12;

            % Create nEditField_5Label
            app.nEditField_5Label = uilabel(app.Defects1to6Panel);
            app.nEditField_5Label.HorizontalAlignment = 'right';
            app.nEditField_5Label.Position = [332 59 25 22];
            app.nEditField_5Label.Text = 'ÿ_n';

            % Create nEditField_5
            app.nEditField_5 = uieditfield(app.Defects1to6Panel, 'numeric');
            app.nEditField_5.Limits = [1e-20 1e-10];
            app.nEditField_5.ValueChangedFcn = createCallbackFcn(app, @dopantFieldValueChanged, true);
            app.nEditField_5.Position = [359 59 53 22];
            app.nEditField_5.Value = 2e-18;

            % Create pEditField_4Label
            app.pEditField_4Label = uilabel(app.Defects1to6Panel);
            app.pEditField_4Label.HorizontalAlignment = 'right';
            app.pEditField_4Label.Position = [332 27 25 22];
            app.pEditField_4Label.Text = 'ÿ_p';

            % Create pEditField_4
            app.pEditField_4 = uieditfield(app.Defects1to6Panel, 'numeric');
            app.pEditField_4.Limits = [1e-20 1e-10];
            app.pEditField_4.ValueChangedFcn = createCallbackFcn(app, @dopantFieldValueChanged, true);
            app.pEditField_4.Position = [359 27 53 22];
            app.pEditField_4.Value = 2e-18;

            % Create DitEditField_9Label
            app.DitEditField_9Label = uilabel(app.Defects1to6Panel);
            app.DitEditField_9Label.HorizontalAlignment = 'right';
            app.DitEditField_9Label.Position = [456 153 18 22];
            app.DitEditField_9Label.Text = 'Dit';

            % Create DitEditField_9
            app.DitEditField_9 = uieditfield(app.Defects1to6Panel, 'numeric');
            app.DitEditField_9.Limits = [10000000 1e+15];
            app.DitEditField_9.ValueChangedFcn = createCallbackFcn(app, @dopantFieldValueChanged, true);
            app.DitEditField_9.Position = [476 153 53 22];
            app.DitEditField_9.Value = 10000000;

            % Create EEtEditField_5Label
            app.EEtEditField_5Label = uilabel(app.Defects1to6Panel);
            app.EEtEditField_5Label.HorizontalAlignment = 'right';
            app.EEtEditField_5Label.Position = [440 123 29 22];
            app.EEtEditField_5Label.Text = 'E-Et';

            % Create EEtEditField_5
            app.EEtEditField_5 = uieditfield(app.Defects1to6Panel, 'numeric');
            app.EEtEditField_5.Limits = [-1 1];
            app.EEtEditField_5.ValueChangedFcn = createCallbackFcn(app, @dopantFieldValueChanged, true);
            app.EEtEditField_5.Position = [476 123 53 22];

            % Create DitEditField_10Label
            app.DitEditField_10Label = uilabel(app.Defects1to6Panel);
            app.DitEditField_10Label.HorizontalAlignment = 'right';
            app.DitEditField_10Label.Position = [441 91 28 22];
            app.DitEditField_10Label.Text = 'ÿDit';

            % Create DitEditField_10
            app.DitEditField_10 = uieditfield(app.Defects1to6Panel, 'numeric');
            app.DitEditField_10.Limits = [-2 2];
            app.DitEditField_10.ValueChangedFcn = createCallbackFcn(app, @dopantFieldValueChanged, true);
            app.DitEditField_10.Position = [476 91 53 22];
            app.DitEditField_10.Value = 1;

            % Create nEditField_6Label
            app.nEditField_6Label = uilabel(app.Defects1to6Panel);
            app.nEditField_6Label.HorizontalAlignment = 'right';
            app.nEditField_6Label.Position = [444 59 25 22];
            app.nEditField_6Label.Text = 'ÿ_n';

            % Create nEditField_6
            app.nEditField_6 = uieditfield(app.Defects1to6Panel, 'numeric');
            app.nEditField_6.Limits = [1e-20 1e-10];
            app.nEditField_6.ValueChangedFcn = createCallbackFcn(app, @dopantFieldValueChanged, true);
            app.nEditField_6.Position = [476 59 53 22];
            app.nEditField_6.Value = 1e-15;

            % Create pEditField_5Label
            app.pEditField_5Label = uilabel(app.Defects1to6Panel);
            app.pEditField_5Label.HorizontalAlignment = 'right';
            app.pEditField_5Label.Position = [444 27 25 22];
            app.pEditField_5Label.Text = 'ÿ_p';

            % Create pEditField_5
            app.pEditField_5 = uieditfield(app.Defects1to6Panel, 'numeric');
            app.pEditField_5.Limits = [1e-20 1e-10];
            app.pEditField_5.ValueChangedFcn = createCallbackFcn(app, @dopantFieldValueChanged, true);
            app.pEditField_5.Position = [476 27 53 22];
            app.pEditField_5.Value = 1e-16;

            % Create typeDropDown_6Label
            app.typeDropDown_6Label = uilabel(app.Defects1to6Panel);
            app.typeDropDown_6Label.HorizontalAlignment = 'right';
            app.typeDropDown_6Label.Visible = 'off';
            app.typeDropDown_6Label.Position = [570 174 29 22];
            app.typeDropDown_6Label.Text = 'type';

            % Create typeDropDown_6
            app.typeDropDown_6 = uidropdown(app.Defects1to6Panel);
            app.typeDropDown_6.Items = {'Don -tail', 'Acc -tail', 'Acc -tophat', 'Don -tophat', 'Acc -gauss', 'Don -gauss'};
            app.typeDropDown_6.ValueChangedFcn = createCallbackFcn(app, @dopantFieldValueChanged, true);
            app.typeDropDown_6.Position = [570 174 99 22];
            app.typeDropDown_6.Value = 'Don -tophat';

            % Create DitEditField_12Label
            app.DitEditField_12Label = uilabel(app.Defects1to6Panel);
            app.DitEditField_12Label.HorizontalAlignment = 'right';
            app.DitEditField_12Label.Position = [563 151 18 22];
            app.DitEditField_12Label.Text = 'Dit';

            % Create DitEditField_12
            app.DitEditField_12 = uieditfield(app.Defects1to6Panel, 'numeric');
            app.DitEditField_12.Limits = [10000000 1e+15];
            app.DitEditField_12.ValueChangedFcn = createCallbackFcn(app, @dopantFieldValueChanged, true);
            app.DitEditField_12.Position = [583 151 53 22];
            app.DitEditField_12.Value = 10000000;

            % Create nEditField_7Label
            app.nEditField_7Label = uilabel(app.Defects1to6Panel);
            app.nEditField_7Label.HorizontalAlignment = 'right';
            app.nEditField_7Label.Position = [551 59 25 22];
            app.nEditField_7Label.Text = 'ÿ_n';

            % Create nEditField_7
            app.nEditField_7 = uieditfield(app.Defects1to6Panel, 'numeric');
            app.nEditField_7.Limits = [1e-20 1e-10];
            app.nEditField_7.ValueChangedFcn = createCallbackFcn(app, @dopantFieldValueChanged, true);
            app.nEditField_7.Position = [583 57 53 22];
            app.nEditField_7.Value = 1e-15;

            % Create pEditField_6Label
            app.pEditField_6Label = uilabel(app.Defects1to6Panel);
            app.pEditField_6Label.HorizontalAlignment = 'right';
            app.pEditField_6Label.Position = [551 27 25 22];
            app.pEditField_6Label.Text = 'ÿ_p';

            % Create pEditField_6
            app.pEditField_6 = uieditfield(app.Defects1to6Panel, 'numeric');
            app.pEditField_6.Limits = [1e-20 1e-10];
            app.pEditField_6.ValueChangedFcn = createCallbackFcn(app, @dopantFieldValueChanged, true);
            app.pEditField_6.Position = [583 25 53 22];
            app.pEditField_6.Value = 1e-16;

            % Create EEtEditField_7Label
            app.EEtEditField_7Label = uilabel(app.Defects1to6Panel);
            app.EEtEditField_7Label.HorizontalAlignment = 'right';
            app.EEtEditField_7Label.Position = [547 123 29 22];
            app.EEtEditField_7Label.Text = 'E-Et';

            % Create EEtEditField_6
            app.EEtEditField_6 = uieditfield(app.Defects1to6Panel, 'numeric');
            app.EEtEditField_6.Limits = [-1 1];
            app.EEtEditField_6.ValueChangedFcn = createCallbackFcn(app, @dopantFieldValueChanged, true);
            app.EEtEditField_6.Position = [583 123 53 22];

            % Create DitEditField_13Label
            app.DitEditField_13Label = uilabel(app.Defects1to6Panel);
            app.DitEditField_13Label.HorizontalAlignment = 'right';
            app.DitEditField_13Label.Position = [547 93 28 22];
            app.DitEditField_13Label.Text = 'ÿDit';

            % Create DitEditField_11
            app.DitEditField_11 = uieditfield(app.Defects1to6Panel, 'numeric');
            app.DitEditField_11.Limits = [-2 2];
            app.DitEditField_11.ValueChangedFcn = createCallbackFcn(app, @dopantFieldValueChanged, true);
            app.DitEditField_11.Position = [582 93 53 22];
            app.DitEditField_11.Value = 1;

            % Create ButtonGroup
            app.ButtonGroup = uibuttongroup(app.Defects1to6Panel);
            app.ButtonGroup.Position = [2 -5 679 30];

            % Create ChargedButton
            app.ChargedButton = uiradiobutton(app.ButtonGroup);
            app.ChargedButton.Text = 'Charged';
            app.ChargedButton.Position = [250 4 68 22];
            app.ChargedButton.Value = true;

            % Create NeutralButton
            app.NeutralButton = uiradiobutton(app.ButtonGroup);
            app.NeutralButton.Text = 'Neutral';
            app.NeutralButton.Position = [317 5 61 22];

            % Create Label_4
            app.Label_4 = uilabel(app.UIFigure);
            app.Label_4.FontSize = 20;
            app.Label_4.Position = [9 606 330 21];
            app.Label_4.Text = '-----------------------------------------';

            % Create Label_5
            app.Label_5 = uilabel(app.UIFigure);
            app.Label_5.FontSize = 20;
            app.Label_5.Position = [9 530 329 24];
            app.Label_5.Text = '-----------------------------------------';

            % Create Label_6
            app.Label_6 = uilabel(app.UIFigure);
            app.Label_6.FontSize = 20;
            app.Label_6.Position = [11 650 330 26];
            app.Label_6.Text = '________________________________';

            % Create Label_7
            app.Label_7 = uilabel(app.UIFigure);
            app.Label_7.HorizontalAlignment = 'right';
            app.Label_7.Position = [982 399 25 22];
            app.Label_7.Text = '#';

            % Create NumSweepEditField_2
            app.NumSweepEditField_2 = uieditfield(app.UIFigure, 'numeric');
            app.NumSweepEditField_2.Limits = [1 50];
            app.NumSweepEditField_2.ValueDisplayFormat = '%.0f';
            app.NumSweepEditField_2.Position = [1012 399 25 22];
            app.NumSweepEditField_2.Value = 20;

            % Create Sn0cmsEditFieldLabel
            app.Sn0cmsEditFieldLabel = uilabel(app.UIFigure);
            app.Sn0cmsEditFieldLabel.HorizontalAlignment = 'right';
            app.Sn0cmsEditFieldLabel.Position = [403 198 63 22];
            app.Sn0cmsEditFieldLabel.Text = 'Sn0 (cm/s)';

            % Create Sn0cmsEditField
            app.Sn0cmsEditField = uieditfield(app.UIFigure, 'numeric');
            app.Sn0cmsEditField.ValueDisplayFormat = '%2.2f';
            app.Sn0cmsEditField.Editable = 'off';
            app.Sn0cmsEditField.Position = [472 198 64 22];

            % Create Sp0cmsEditFieldLabel
            app.Sp0cmsEditFieldLabel = uilabel(app.UIFigure);
            app.Sp0cmsEditFieldLabel.HorizontalAlignment = 'right';
            app.Sp0cmsEditFieldLabel.Position = [404 167 63 22];
            app.Sp0cmsEditFieldLabel.Text = 'Sp0 (cm/s)';

            % Create Sp0cmsEditField
            app.Sp0cmsEditField = uieditfield(app.UIFigure, 'numeric');
            app.Sp0cmsEditField.ValueDisplayFormat = '%2.2f';
            app.Sp0cmsEditField.Editable = 'off';
            app.Sp0cmsEditField.Position = [473 167 64 22];

            % Create J0sfAcm2EditFieldLabel
            app.J0sfAcm2EditFieldLabel = uilabel(app.UIFigure);
            app.J0sfAcm2EditFieldLabel.HorizontalAlignment = 'right';
            app.J0sfAcm2EditFieldLabel.Position = [393 137 73 22];
            app.J0sfAcm2EditFieldLabel.Text = 'J0s (fA/cm2)';

            % Create J0sfAcm2EditField
            app.J0sfAcm2EditField = uieditfield(app.UIFigure, 'numeric');
            app.J0sfAcm2EditField.ValueDisplayFormat = '%2.2f';
            app.J0sfAcm2EditField.Editable = 'off';
            app.J0sfAcm2EditField.Position = [472 137 64 22];

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = srv_app_matlab_raw

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end