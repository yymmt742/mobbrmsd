character(*), parameter  :: title_default       = 'element_set default'
character(*), parameter  :: description_default = 'default element set.'
type(element), parameter :: element_set_default(10) = [&
    element('H',  'H',  'H',   1, 1.00797_RK, 0.0_RK),&
    element('He', 'He', 'He',  2, 4.00260_RK, 0.0_RK),&
    element('Li', 'Li', 'Li',  3, 6.94100_RK, 0.0_RK),&
    element('Be', 'Be', 'He',  4, 9.01218_RK, 0.0_RK),&
    element('B',  'B',  'H',   5, 10.0110_RK, 0.0_RK),&
    element('C',  'C',  'H',   6, 12.0110_RK, 0.0_RK),&
    element('N',  'N',  'H',   7, 14.0067_RK, 0.0_RK),&
    element('O',  'O',  'H',   8, 15.9994_RK, 0.0_RK),&
    element('F',  'F',  'H',   9, 18.9984_RK, 0.0_RK),&
    element('Ne', 'Ne', 'He', 10, 20.1790_RK, 0.0_RK)&
    ]
