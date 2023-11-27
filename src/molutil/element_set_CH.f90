character(*), parameter  :: title_CH       = 'element_set HC'
character(*), parameter  :: description_CH = 'Simple element set by Hydrogen and Carbon.'
type(element), parameter :: element_set_CH(2) = [&
    element('H',  'H',  'H',   1, 1.00797_RK, 0.0_RK),&
    element('C',  'C',  'H',   6, 12.0110_RK, 0.0_RK) &
    ]
