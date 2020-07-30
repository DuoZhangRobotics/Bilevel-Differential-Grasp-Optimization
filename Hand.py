import numpy as np
import vtk
from hand import Hand, vtk_render, vtk_add_from_hand, Link

if __name__ == '__main__':
    path = 'hand/ShadowHand/'
    scale = 0.01
    use_eigen = True
    hand = Hand(path, scale, use_joint_limit=True, use_quat=True, use_eigen=use_eigen)
    # if hand.use_eigen:
    #     dofs = 1000 * np.ones(hand.eg_num)
    # else:
    #     dofs = 1000 * np.ones(hand.nr_dof())
    # hand.forward_kinematics(1000 * np.ones(hand.extrinsic_size), dofs)
    # hand.value_check(10)
    # hand.grad_check(2)
    # # hand.write_limits()
    # renderer = vtk.vtkRenderer()
    # vtk_add_from_hand(hand, renderer, scale)
    # vtk_render(renderer, axes=False)
