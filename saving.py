# draw for link
def draw(self, save=False, path=None, idx=None):
    if self.mesh:
        ret = copy.deepcopy(self.mesh).apply_transform(self.joint_transform)
        if save:
            obj = trimesh.exchange.obj.export_obj(ret)
            with open(os.path.join(path, str(idx) + '.obj'), 'w') as f:
                f.write(obj)
    else:
        ret = None
    idx += 1
    if self.children:
        for c in self.children:
            if ret:
                tmp_ret, idx = c.draw(save, path, idx)
                ret += tmp_ret
            else:
                ret, idx = c.draw(save, path, idx)
    return ret, idx


# draw for hand
def draw(self, scale_factor=1, show_to_screen=True, save=False, path=None):
    if save and not os.path.exists(path):
        os.makedirs(path)

    mesh, _ = self.palm.draw(save, path, 0)
    mesh.apply_scale(scale_factor)
    if show_to_screen:
        mesh.show()
    return mesh


hand.draw(scale_factor=1, show_to_screen=False, save=True, path='transformed_objs')