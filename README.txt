HW1 
Name: Jing Qiu
Penn ID: 43429978

Hi! This is the introduction document for the Jello Simulation project. First I will answer the five questions on the homework website.
1.The effect of Ks and Kd. Well, in plain words, a larger Ks gives us a stiffer Jello, which is harder to deform. Kd describes the damping effect. A bigger Kd will decrement the response time for the Jello to reach a stable state.From the mathmetical point of view,the internal forces (spring forces) between particles are calculated as the sum of spring elastic forces and damping forces. Ks is the stiffness contant for the elastic force which is given by the hook's law, and the Kd is the damping ratio of the damping force.
2. Drawbacks: Requires impractically small time steps. When the time step increment, the particles are more likely to reach an invalid position inside the obstacle.
   Benefits: Easy to implement. The virtual spring force is a pretty good simulation of real collision external forces.
3. Mass-Spring System
4. In plain words, Explicit Euler methods calculate the state of a system at a later time from the state of the system at the current time, while implicit Euler methods find a solution by solving an equation involving both the current state of the system and the later one. 
5. Yes, my Jello does! Check out the video output.

For basic feature, using the default parameters with StructuralKs = 5000, shearKs = 4000,bendKs = 4000,StructuralKd = 10, shearKd = 10,bendKd = 10. Make sure that in the JelloMesh::Update function,the ResolveDragging function is commented out, since this is for the extra credit.

In the JelloMesh.cpp,you can change the default integration type in the first funtion JelloMesh::JelloMesh(), by changing m_integrationType(JelloMesh::RK4) to (JelloMesh::EULER) or JelloMesh::MIDPOINT . Since I didn't add many additional springs to the Jello, I manually make the Jello "softer" for EULER and MIDPOINT simulation by changing the initialize parameters. Otherwise, the Jello is likely to explode under the default settings. To change the initial parameters, just commented out the default one and uncomment the StructuralKs = 1000, shearKs = 1000,bendKs = 1000,StructuralKd = 5, shearKd = 5,bendKd = 5.
(Note I set the VirtualKd and VirtualKs as global variables used in the contact and collision senarios instead of using penalyKs and penaltyKd, so don't get panic when seeing penaltyKs or penaltyKd is zero)

The environment or in other word, the world set up can be changed in the main.cpp file, World theWorld("worlds/ground.xml") can be changed to World theWorld("worlds/cylinders.xml") for the second set up. As for extra credit 1, whereas the Jello interact with sphere and cube, use World theWorld("worlds/cube_sphere.xml")

Ok, so far the implementation is kind of easy and basic. I'll talk a little bit more about the extra credit 2 and 3, which are the user dragging particles and user generate force. First of all, I wanna thank the TA Sijie, who helped me a lot on implementing this feature.
The feature still has a lot of problems unsolved, but at least I can provide a descent video with a user external force dragging the Jello in the output folder. The main problem is the stability of the Jello. Hopefully, this can be solved in the future by adding more and more springs to the Jello or implementing the implicit intergration method, which I really don't have time for that.

So, let's get started. Firstly, we need to make the Jello stiffer. Change the initialize parameters to StructuralKs = 50000, shearKs = 40000,bendKs = 40000,StructuralKd = 1000, shearKd = 1000,bendKd = 1000. The next thing to do is uncomment the ResolveDragging function in the JelloMesh::Update (This is important!!!) Without this step, you won't see anything happen when you click you mouse. And make sure the integration type is RK4, which is the best for now.
Now we are all set and ready to simulate. Run the program, and keep an eye on the console window as well, since I output some messages there for you to see. If you click on the Jello cube, the console window will show up messages: "Got it! Pick up particle xx", xx is the index of the particle, this tells you which particle you have selected by the mouse. Next, if you still click on the Jello, you will probably change the selected particle but no additional force is genreated. Now, if you press on the outside of the Jello, the console window will print out two set of coordinates. First set is the position of the chosen particle and the second set is the desired position specified by the mouse position outside the Jello, then you will probably see the Jello being dragged. Again, because of the stability problem. It's likely to explode. 

I won't go too deep on how the dragging is realized. The main point is the function vec3 GetOGLPos(int x, int y) in main.cpp which takes on the mouse position x,y and translates them into a world position (PosX1,PosY1,PosZ1) and return it. I will give credit to the following website for the implementation of this function.  http://nehe.gamedev.net/article/using_gluunproject/16013/

So, simply speaking, I pass the mouse coordinates in the world from main.cpp to JelloMesh.cpp, and compare them with every particle in JelloMesh, and return the index of the selected particle to the main.cpp then we pass relevant parameters to JelloMesh.cpp calling the dragJello function to register dragging action. And finally in the ResolveDragging function we solve the dragging force. 
(Notice that: in order to keep track of the user generated force, I create a another field in the Particle class called userforce, used to store the dragging force for each particle)



Thank you! And enjoy!

Reference: TA Sijie and Issam
Website: http://nehe.gamedev.net/article/using_gluunproject/16013/
         http://www.songho.ca/opengl/gl_transform.html
         http://www.opengl.org/resources/libraries/glut/spec3/node51.html
Books on Linear Algebra         
Classmate: Naichen and Zimeng   

 



