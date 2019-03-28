from coopr.pyomo import *


def phase_1_2_integrate(phase1_inst, phase2_inst):
    
    n_tours = sum(phase1_inst.TourType[i,t].value for (i,t) in phase1_inst.okTourType)
        
    phase2_inst.n_tours = Param(initialize=n_tours)
        
    phase2_inst.TOURS = RangeSet(1,phase2_inst.n_tours)
      
    for (i,j,w,k,t) in phase1_inst.okShifts:
        phase2_inst.Shift[i,j,w,k,t].value = phase1_inst.Shift[i,j,w,k,t].value

    for (i,t) in phase1_inst.TourType_index:
        phase2_inst.TourType[i,t].value = phase1_inst.TourType[i,t].value
                 
    for (i,t,j,w) in phase1_inst.DailyTourType_index:
        phase2_inst.DailyTourType[i,t,j,w].value = phase1_inst.DailyTourType[i,t,j,w].value
                                                                         
    for (i,t,k,j,w) in phase1_inst.ok_daily_shift_index:
        phase2_inst.DailyShiftWorked[i,t,k,j,w].value = phase1_inst.DailyShiftWorked[i,t,k,j,w].value
                                                                         
    for (i,t,p) in phase1_inst.ok_weekenddaysworked_index:
        phase2_inst.WeekendDaysWorked[i,t,p].value = phase1_inst.WeekendDaysWorked[i,t,p].value                                                                 
      
        # Initialize the phase 2 instance
        
    print n_tours
