% Various examples of torque free motion

binet_construction
%%
torque_free_rigid_body('w0', [2,0.5,0.5])

%%
torque_free_rigid_body('w0', [0.5,0.5,2])

%%
torque_free_rigid_body('w0', [0,2,0]+1e-3)

%%
torque_free_rigid_body('poinsot', true)