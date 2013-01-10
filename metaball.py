import bpy
import os

x=0.8
y=1
z=1.2
i=1

mb=bpy.data.metaballs.new(name="mymb"+str(i))
el=bpy.data.metaballs['mymb'+str(i)].elements.new('ELLIPSOID')
ob=bpy.data.objects.new('myel0'+str(i),mb)
bpy.context.scene.objects.link(ob)
el.size_x=x
el.size_y=y
el.size_z=z
bpy.data.metaballs['mymb'+str(i)].resolution=0.01
bpy.context.scene.objects.active=ob
ob.select=True
bpy.ops.object.convert(target='MESH', keep_original=False)
bpy.ops.object.editmode_toggle()
bpy.ops.mesh.select_all(action='TOGGLE')
bpy.ops.mesh.quads_convert_to_tris()
bpy.ops.object.modifier_add(type='DECIMATE')
bpy.ops.object.editmode_toggle()
bpy.data.objects['myel0'+str(i)+'.001'].modifiers['Decimate'].ratio=0.5
bpy.ops.object.modifier_apply(apply_as='DATA',modifier="Decimate")

file = open("ellipsoid("+str(x)+","+str(y)+","+str(z)+").txt", "w", encoding="utf8", newline="\n")
mesh = bpy.context.object.data
file.write('%d %d\n' % (len(mesh.vertices),len(mesh.tessfaces)))
for k in range(len(mesh.vertices)):
	co=mesh.vertices[k].co
	file.write('%.6f %.6f %.6f\n' % (co[0],co[1],co[2]))
for k in range(len(mesh.tessfaces)):
	fa=mesh.tessfaces[k].vertices
	file.write('%d %d %d\n' % (fa[0],fa[1],fa[2]))
file.close()
