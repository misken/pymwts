fin = open("simple_w0_phase2_results.yml","r")
sol = yaml.load(fin)

# Get name of the problem
sol['Problem'][0]['Name']

# Get solution status
sol['Solution'][1]['Status']

# Get variable information
sol['Solution'][1]['Variable']['TourWkendDof(1,59,13,1)']['Id']
sol['Solution'][1]['Variable']['TourWkendDof(1,59,13,1)']['Value']
