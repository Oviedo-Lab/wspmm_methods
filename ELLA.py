
from ELLA.ELLA import model_beta, model_null, loss_ll, ELLA

ella_demo = ELLA(dataset='demo1', adam_learning_rate_min=1e-2, max_iter=1000)

ella_demo.load_data(data_path='input/mini_demo_data.pkl')
ella_demo.register_cells()
ella_demo.nhpp_prepare()

