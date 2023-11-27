character(*), parameter  :: title_HEAVY       = 'element_set HEAVY'
character(*), parameter  :: description_HEAVY = 'element set of HEAVY atoms.'
type(element), parameter :: element_set_HEAVY(8) = [&
    element('Li', 'Li', 'Li',  3, 6.94100_RK, 0.0_RK),&
    element('Be', 'Be', 'He',  4, 9.01218_RK, 0.0_RK),&
    element('B',  'B',  'H',   5, 10.0110_RK, 0.0_RK),&
    element('C',  'C',  'H',   6, 12.0110_RK, 0.0_RK),&
    element('N',  'N',  'H',   7, 14.0067_RK, 0.0_RK),&
    element('O',  'O',  'H',   8, 15.9994_RK, 0.0_RK),&
    element('F',  'F',  'H',   9, 18.9984_RK, 0.0_RK),&
    element('Ne', 'Ne', 'He', 10, 20.1790_RK, 0.0_RK)&
    ]
