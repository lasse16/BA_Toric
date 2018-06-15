/*
 * Created by SharpDevelop.
 * User: Lasse Haffke
 * Date: 14.06.2018
 * Time: 11:13
 * 
 *
*  To change this template use Tools | Options | Coding | Edit Standard Headers.
 */
using System;
using UnityEngine;



/// <summary>
/// A representtaation of an angle, always under 360 degrees and with a easy toRad function.
/// </summary>
public class FixAngle
{
	private float _angle;
	
	public FixAngle(float angle)
	{
		_angle = angle % 360;
	}
	
	public float angle()
	{
		return _angle;
	}
	
	public float toRad(){
		return _angle * Deg2Rad;
	}
}

